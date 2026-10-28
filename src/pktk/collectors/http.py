"""Small public-file collector: robots, retries, atomic cache, digest and audit log."""
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit
from urllib.robotparser import RobotFileParser
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

AGENT = 'pktk-research/0.1'

def now():
    return datetime.now(timezone.utc).isoformat()

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def log(path, **event):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as f:
        f.write(json.dumps({'time': now(), **event}) + '\n')

class Collector:
    def __init__(self, root, offline=False):
        self.root = Path(root)
        self.raw = self.root / 'data/raw'
        self.raw.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.session = requests.Session()
        self.session.headers['User-Agent'] = AGENT
        retry = Retry(total=4, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504], respect_retry_after_header=True)
        self.session.mount('https://', HTTPAdapter(max_retries=retry))
        self.robots = {}

    def allowed(self, url):
        u = urlsplit(url)
        origin = f'{u.scheme}://{u.netloc}'
        if origin not in self.robots:
            r = self.session.get(origin + '/robots.txt', timeout=60)
            rp = RobotFileParser()
            if r.status_code in (404, 410):
                rp.parse([])
            else:
                r.raise_for_status()
                rp.parse(r.text.splitlines())
            self.robots[origin] = rp
            log(self.raw/'collection.jsonl', action='robots', url=origin+'/robots.txt', status=r.status_code,
                sha256=hashlib.sha256(r.content).hexdigest())
        rp = self.robots[origin]
        if not rp.can_fetch(AGENT, url):
            raise PermissionError(f'robots.txt disallows {url}')
        time.sleep(max(0.3, rp.crawl_delay(AGENT) or 0))

    def fetch(self, name, url, version, expected_sha256=None, json_body=None, allowed_statuses=()):
        target = self.raw/name
        manifest = target.with_suffix(target.suffix+'.json')
        try:
            if target.exists() and manifest.exists():
                meta = json.loads(manifest.read_text())
                digest = sha256(target)
                if meta.get('http_status',200)>=400 and meta['http_status'] not in allowed_statuses:
                    raise ValueError('Cached HTTP error not allowed for this request')
                if meta['url'] != url or meta['version'] != version or digest != meta['sha256'] or meta.get('request_body') != json_body:
                    raise ValueError(f'Cache provenance or integrity mismatch: {name}')
                if expected_sha256 and digest != expected_sha256:
                    raise ValueError(f'Pinned checksum mismatch: {name}')
                log(self.raw/'collection.jsonl', action='cache_hit', file=name)
                return target, meta
            if self.offline:
                raise FileNotFoundError(f'No verified cache for {name}; run online first')
            self.allowed(url)
            # No authentication or cross-host redirects are followed implicitly.
            if json_body is None:
                r = self.session.get(url, timeout=120, allow_redirects=False)
            else:
                r = self.session.post(url, json=json_body, timeout=120, allow_redirects=False)
            if r.status_code not in allowed_statuses:
                r.raise_for_status()
            if r.is_redirect:
                raise PermissionError('Redirect requires endpoint review: '+r.headers.get('Location',''))
            digest = hashlib.sha256(r.content).hexdigest()
            if expected_sha256 and digest != expected_sha256:
                raise ValueError(f'Pinned checksum mismatch: {name}')
            temp = target.with_suffix(target.suffix+'.part')
            temp.write_bytes(r.content)
            temp.replace(target)
            meta = dict(url=url, version=version, retrieved_at=now(), sha256=digest,
                        bytes=len(r.content), http_status=r.status_code, etag=r.headers.get('ETag'), last_modified=r.headers.get('Last-Modified'))
            if json_body is not None:
                meta.update(method='POST',request_body=json_body)
            manifest.write_text(json.dumps(meta, indent=2)+'\n')
            log(self.raw/'collection.jsonl', action='download', file=name, **meta)
            return target, meta
        except Exception as e:
            log(self.raw/'collection.jsonl', action='failed', url=url, error=str(e))
            raise

    def pages(self, url):
        """For JSON APIs exposing results/next; each page is cached, loops rejected."""
        seen = set()
        while url:
            if url in seen:
                raise ValueError('Pagination loop')
            seen.add(url)
            path, _ = self.fetch(hashlib.sha256(url.encode()).hexdigest()+'.json', url, 'api-snapshot')
            page = json.loads(path.read_text())
            yield page['results']
            url = page.get('next')
