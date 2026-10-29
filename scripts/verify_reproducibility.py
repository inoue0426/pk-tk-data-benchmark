"""Replay the complete cached pipeline and verify artifacts against the preceding run."""
import hashlib,json,os,platform,re,subprocess,sys
from datetime import datetime, timezone
from pathlib import Path
root=Path(__file__).resolve().parents[1]
env=dict(os.environ,PYTHONPATH=str(root/'src'))
report_names=['summary.json','expansion_summary.json','composite_audit.csv','composite_audit.json',
              'cebs_row_selection.csv','cebs_expansion_manifest.json','expansion_overlap.json',
              'stage_snapshots.json','qc_counts.json','final_report.md','benchmark_feasibility.md',
              'scientific_analysis.md','one_page_summary.md','research_manifest.json']
files=sorted(list((root/'data/processed').glob('*.csv'))+list((root/'data/processed').glob('*.json'))+
             [root/'reports'/name for name in report_names]+list((root/'results/research').glob('*.csv'))+
             list((root/'results/research').glob('*.json'))+list((root/'results/figures').glob('*.png'))+
             list((root/'results/figures').glob('*.pdf'))+list((root/'results/figures').glob('*.svg')))
def hashes():return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
before=hashes()
subprocess.run([sys.executable,'-m','pktk.pipeline','--offline'],cwd=root,env=env,check=True)
after=hashes();differences=[p for p in before if before[p]!=after[p]]
tests=subprocess.run([sys.executable,'-m','pytest','-q'],cwd=root,env=env,capture_output=True,text=True)
print(tests.stdout);print(tests.stderr)
subprocess.run([sys.executable,'-m','compileall','-q','src','scripts'],cwd=root,check=True)
subprocess.run(['git','diff','--check'],cwd=root,check=True)
passed=re.search(r'(\d+) passed',tests.stdout)
validation=dict(validation_date=datetime.now(timezone.utc).date().isoformat(),offline_full_pipeline='passed',
    byte_identical_data_and_report_artifacts=len(files)-len(differences),differences=differences,
    tests_passed=int(passed[1]) if passed else 0,tests_exit_code=tests.returncode,
    python=platform.python_version(),platform=platform.platform(),git_diff_check='passed',compileall='passed',
    warnings='openpyxl warns about unsupported extensions. Source workbooks are read only and checksum verified.',
    model_scope='Expanded empirical scaling, power-law, disjoint pooled/shrinkage and identifiable nonlinear models; no complex molecular ML',artifact_sha256=after)
(root/'reports/validation.json').write_text(json.dumps(validation,indent=2)+'\n')
(root/'reports/validation.md').write_text('# Reproducibility validation\n\n'+
    f"Complete offline collection → harmonization → QC replay executed. {len(files)-len(differences)}/{len(files)} selected data/report artifacts are byte-identical to the preceding run. Differences: {differences}.\n\n"+
    f"Tests: {validation['tests_passed']} passed; pytest exit code {tests.returncode}. Compilation and git diff whitespace checks passed. See validation.json for hashes. Expanded empirical benchmarks and all figures were regenerated; no complex molecular ML was trained.\n\n"+
    'The replay uses locally cached source bytes and identity API snapshots. It does not prove future endpoint availability or bitwise equivalence of fresh external API responses. Run `PYTHONPATH=src python scripts/verify_reproducibility.py` after an online full/expansion run.\n')
if differences or tests.returncode:raise SystemExit('Reproducibility or tests failed; see reports/validation.json')
