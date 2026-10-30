# Dataset A reproducibility validation

Full Dataset A offline collection → identity → harmonization → QC → feasibility/report replay executed. 14/14 A artifacts are byte-identical; 44/44 B data/results/figure artifacts are unchanged. No strict baseline rerun.

Pytest exit code: 0.

```text
........................................................                 [100%]
=============================== warnings summary ===============================
tests/test_cebs_reviewed.py: 9 warnings
tests/test_expansion.py: 12 warnings
  /Users/inouey2/code/pk-tk-data-benchmark/.venv/lib/python3.14/site-packages/openpyxl/worksheet/_reader.py:329: UserWarning: Unknown extension is not supported and will be removed
    warn(msg)

tests/test_molecular_ml.py: 62 warnings
  /Users/inouey2/code/pk-tk-data-benchmark/.venv/lib/python3.14/site-packages/torch_geometric/inspector.py:433: DeprecationWarning: Failing to pass a value to the 'type_params' parameter of 'typing._eval_type' is deprecated, as it leads to incorrect behaviour when calling typing._eval_type on a stringified annotation that references a PEP 695 type parameter. It will be disallowed in Python 3.15.
    return typing._eval_type(value, _globals, None)  # type: ignore

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
56 passed, 83 warnings in 14.76s
```

Hashes and limitations: molecular_validation.json. Cached replay verifies stored bytes, not future API availability.
