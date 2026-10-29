# Dataset A reproducibility validation

Full Dataset A offline collection → identity → harmonization → QC → feasibility/report replay executed. 14/14 A artifacts are byte-identical; 44/44 B data/results/figure artifacts are unchanged. No strict baseline rerun.

Pytest exit code: 0.

```text
..........................................                               [100%]
=============================== warnings summary ===============================
tests/test_cebs_reviewed.py: 9 warnings
tests/test_expansion.py: 12 warnings
  /Users/inouey2/code/pk-tk-data-benchmark/.venv/lib/python3.14/site-packages/openpyxl/worksheet/_reader.py:329: UserWarning: Unknown extension is not supported and will be removed
    warn(msg)

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
42 passed, 21 warnings in 8.53s
```

Hashes and limitations: molecular_validation.json. Cached replay verifies stored bytes, not future API availability.
