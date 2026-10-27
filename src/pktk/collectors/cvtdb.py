import tarfile
import warnings
import rdata

URL = 'https://cran.r-project.org/src/contrib/invivoPKfit_2.0.2.tar.gz'

def collect(collector, checksum):
    path, meta = collector.fetch('invivoPKfit_2.0.2.tar.gz', URL, 'invivoPKfit-2.0.2/CvTdb-2025-08-12', checksum)
    out = collector.root/'data/interim'
    out.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path) as archive:
        # Read only named members; never unpack an untrusted archive onto disk.
        member = out/'cvtdb_original.rda'
        member.write_bytes(archive.extractfile('invivoPKfit/data/cvtdb_original.rda').read())
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='Missing constructor for R class')
            data = rdata.read_rda(member)['cvtdb_original']
        (out/'invivoPKfit_DESCRIPTION.txt').write_bytes(archive.extractfile('invivoPKfit/DESCRIPTION').read())
    data.to_csv(out/'cvtdb_original.csv', index=False)
    return data, meta
