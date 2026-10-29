"""Replay Dataset A without rerunning or changing Dataset B numerical benchmarks."""
import hashlib,json,os,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
ENV=dict(os.environ,PYTHONPATH=str(ROOT/'src'))

def hashes(paths):
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}

def main():
    a=list((ROOT/'data/processed/molecular').glob('*'))+[ROOT/'reports'/n for n in ['molecular_summary.json','molecular_manifest.json','molecular_collection_manifest.json','molecular_pk_expansion.md','molecular_qc.md','molecular_description_audit.csv']]
    if not (ROOT/'reports/molecular_summary.json').exists():raise SystemExit('Run pktk.molecular online first')
    b=[ROOT/'data/processed'/n for n in ['high_confidence_endpoints.csv','crossdose_series_expanded.csv']]+list((ROOT/'results/research').glob('*'))+list((ROOT/'results/figures').glob('*'))+[ROOT/'reports/summary.json']
    before_a=hashes(a);before_b=hashes(b)
    previous=json.loads((ROOT/'reports/validation.json').read_text())['artifact_sha256']
    historical_differences=[p for p,h in before_b.items() if p in previous and previous[p]!=h]
    if historical_differences:raise SystemExit('Dataset B differs from previous verified snapshot: '+str(historical_differences))
    subprocess.run([sys.executable,'-m','pktk.molecular','--offline'],cwd=ROOT,env=ENV,check=True)
    after_a=hashes(a);after_b=hashes(b)
    differences_a=[p for p in before_a if before_a[p]!=after_a[p]]
    differences_b=[p for p in before_b if before_b[p]!=after_b[p]]
    tests=subprocess.run([sys.executable,'-m','pytest','-q'],cwd=ROOT,env=ENV,capture_output=True,text=True)
    print(tests.stdout);print(tests.stderr)
    subprocess.run([sys.executable,'-m','compileall','-q','src','scripts'],cwd=ROOT,check=True)
    subprocess.run(['git','diff','--check'],cwd=ROOT,check=True)
    result=dict(date=datetime.now(timezone.utc).isoformat(),dataset_a_byte_identical=len(a)-len(differences_a),dataset_a_artifacts=len(a),
      dataset_a_differences=differences_a,dataset_b_unchanged_artifacts=len(b)-len(differences_b),dataset_b_differences=differences_b,
      pytest_exit_code=tests.returncode,pytest_output=tests.stdout,artifact_sha256=after_a,dataset_b_preservation_sha256=after_b,
      strict_baselines_rerun=False,dataset_b_matches_previous_snapshot=True,compileall='passed',git_diff_check='passed')
    (ROOT/'reports/molecular_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'reports/molecular_validation.md').write_text(f'# Dataset A reproducibility validation\n\nFull Dataset A offline collection → identity → harmonization → QC → feasibility/report replay executed. {result["dataset_a_byte_identical"]}/{len(a)} A artifacts are byte-identical; {result["dataset_b_unchanged_artifacts"]}/{len(b)} B data/results/figure artifacts are unchanged. No strict baseline rerun.\n\nPytest exit code: {tests.returncode}.\n\n```text\n'+tests.stdout+'```\n\nHashes and limitations: molecular_validation.json. Cached replay verifies stored bytes, not future API availability.\n')
    if differences_a or differences_b or tests.returncode:raise SystemExit('Verification failed; inspect molecular_validation.json')

if __name__=='__main__':main()
