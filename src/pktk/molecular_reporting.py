"""Separate verified Dataset A coverage from original-protocol Dataset B counts."""
import json
from pathlib import Path
import pandas as pd
from pktk.collectors.http import sha256

def table(rows,cols):
    return '\n'.join(['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']+
      ['| '+' | '.join(str(r.get(c,'')) for c in cols)+' |' for r in rows])

def reports(root,x,good,metadata,manifest,pmeta,feasibility,candidates,conflicts):
    from pktk.molecular import dump
    audit_spec=json.loads((root/'configs/molecular/description_audit.json').read_text())
    expected={r['source_record_id']:r['expected_eligible_A'] for r in audit_spec['records']}
    audit=x[x.source_record_id.isin(expected)].copy()
    if len(audit)!=len(expected) or any(bool(r.eligible_A)!=expected[r.source_record_id] for _,r in audit.iterrows()):
        raise ValueError('Reviewed source metadata decisions changed; audit before proceeding')
    cols=['source_record_id','assay_id','publication_id','endpoint','dose_mg_kg','matrix','route','eligible_A','exclusion_reasons','description']
    audit=audit[cols];audit['review_level']=audit_spec['review_level'];audit['license']='CC-BY-SA-3.0'
    audit.to_csv(root/'reports/molecular_description_audit.csv',index=False)
    summary_b=json.loads((root/'reports/summary.json').read_text())
    auc=set(good[good.endpoint=='auc'].inchikey);cmax=set(good[good.endpoint=='cmax'].inchikey)
    strict=good[good.source=='strict_b']
    summary_a=dict(unique_structured_compounds=good.inchikey.nunique(),unique_connectivity_groups=good.connectivity_key.nunique(),
      compounds_with_auc=len(auc),compounds_with_cmax=len(cmax),compounds_with_both_auc_cmax=len(auc&cmax),
      compounds_with_clearance=good[good.endpoint=='clearance'].inchikey.nunique(),
      compounds_with_volume=good[good.endpoint=='volume_distribution'].inchikey.nunique(),
      compounds_with_half_life=good[good.endpoint=='half_life'].inchikey.nunique(),
      collected_endpoint_records=len(x),eligible_endpoint_records=len(good),excluded_endpoint_records=int((~x.eligible_A).sum()),
      known_source_publications=good[good.publication_id.ne('')].publication_id.nunique(),
      independent_experimental_studies=None,original_protocol_verified_studies=strict.original_study_id.nunique(),
      compounds_with_verified_comparable_doses=strict.inchikey.nunique(),
      dose_normalized_records=int(good.dose_normalized.eq('yes').sum()),
      explicit_dose_records=int(good.dose_mg_kg.ne('').sum()),
      target_500_obtained=bool(good.inchikey.nunique()>=500))
    b=dict(validated_series=summary_b['primary_comparable_series'],compounds=summary_b['primary_comparable_compound_identifiers'],
      compounds_with_structures=summary_b['primary_compounds_with_resolved_structures'],
      independent_studies=summary_b['primary_independent_studies'],compounds_two_comparable_doses=summary_b['compounds_at_least_two_comparable_doses'],
      compounds_three_comparable_doses=summary_b['compounds_at_least_three_comparable_doses'],
      auc_series=summary_b['primary_auc_eligible_series'],cmax_series=summary_b['primary_cmax_eligible_series'],
      both_endpoint_series=summary_b['primary_both_endpoints_eligible_series'],new_validated_series=0,
      studies_with_reliable_original_dose_records=summary_b['primary_independent_studies'])
    # B counts include unresolved structures; A's structured intersection is explicitly smaller.
    from pktk.benchmark.episodes import endpoint_table
    e=pd.read_csv(root/'data/processed/high_confidence_endpoints.csv',keep_default_na=False)
    s=pd.read_csv(root/'data/processed/crossdose_series_expanded.csv',keep_default_na=False)
    t=endpoint_table(e,s)
    ba=set(t[t.endpoint.eq('auc')].molecular_group_id);bc=set(t[t.endpoint.eq('cmax')].molecular_group_id)
    b.update(compounds_with_auc=len(ba),compounds_with_cmax=len(bc),compounds_with_both_auc_cmax=len(ba&bc))
    three_ids=set(s[s.primary_eligible & s.n_doses.ge(3)].molecular_group_id) if 'n_doses' in s else set()
    # Use original identity mapping to avoid implying all Dataset A endpoints form dose series.
    structured_ids=set(e[e.canonical_smiles.ne('')].molecular_group_id)
    summary_a['compounds_with_three_verified_comparable_doses']=len(three_ids&structured_ids) if three_ids else None
    sources=[]
    for source,g in x.groupby('source',sort=True):
        q=g[g.eligible_A];sources.append(dict(source=source,collected_records=len(g),eligible_records=len(q),compounds=q.inchikey.nunique(),
           auc_compounds=q[q.endpoint.eq('auc')].inchikey.nunique(),cmax_compounds=q[q.endpoint.eq('cmax')].inchikey.nunique(),
           source_publications=q[q.publication_id.ne('')].publication_id.nunique()))
    reasons=x.loc[~x.eligible_A,'exclusion_reasons'].str.split(';').explode().value_counts().to_dict();reasons.pop('',None)
    summary=dict(dataset_a=summary_a,dataset_b=b,source_coverage=sources,exclusion_flags=reasons,
      identity_conflict_identifiers=len(conflicts),unverified_dose_candidates=len(candidates),
      limitations=['Source publications are not counts of independently verified experiments.',
       'AUC/Cmax coverage includes differently defined and dose-normalized values; these are never pooled.',
       'Source endpoint aggregates and original concentration/time observations remain in separate files.',
       'All raw/derived data remain local; ChEMBL CC BY-SA and PKSmart CC BY obligations are preserved.'])
    dump(root/'reports/molecular_summary.json',summary)
    manifests=manifest['downloads']+[dict(file='pksmart_supplement.xlsx',**pmeta)]
    out=root/'data/processed/molecular'
    dump(root/'reports/molecular_manifest.json',dict(downloads=manifests,
      inputs={str(p.relative_to(root)):sha256(p) for p in [root/'configs/molecular/sources.json',root/'configs/molecular/description_audit.json',root/'data/processed/high_confidence_endpoints.csv',root/'data/processed/crossdose_series_expanded.csv']},
      outputs={str(p.relative_to(root)):sha256(p) for p in sorted(out.iterdir()) if p.suffix in ['.csv','.json']},
      chembl_query_total=manifest['query_total_records'],chembl_collected=manifest['collected_records'],bounded_snapshot=manifest['truncated']))
    report=f'''# Dataset A molecular PK expansion — 2026-10-29

## Verified result and dataset boundary

**Dataset A contains {summary_a['unique_structured_compounds']:,} unique standardized full-InChIKey structures with numerical preclinical PK endpoints** ({summary_a['unique_connectivity_groups']:,} connectivity groups). The 500-structure target is **{'met' if summary_a['target_500_obtained'] else 'not met'}**. These are source-curated numerical measurements with explicit in vivo evidence and RDKit-validated structures, not an assertion that every original animal protocol has been reviewed.

Dataset B remains **{b['validated_series']} validated series, {b['compounds']} compounds and {b['independent_studies']} independently reviewed source studies**. No new B series were promoted and its numerical benchmark was not rerun. {b['compounds_with_structures']} B compounds have resolved structures and contribute to A. Records missing structures are retained in the exclusion ledger. A compound appearing in A never establishes cross-dose eligibility.

## Coverage

{table([dict(metric=k,count=v) for k,v in summary_a.items()],['metric','count'])}

`independent_experimental_studies = null` means unknown, not zero. Dataset A has {summary_a['known_source_publications']} deduplicated publication/source-study proxies with eligible endpoints; compilation rows lacking original citations are not counted as additional studies. Only the inherited {b['independent_studies']} Dataset B studies have original-protocol validation. A numeric dose extracted from an abstracted assay description is not a reliable original protocol.

Dataset B: {b['compounds_two_comparable_doses']} compounds have ≥2 comparable doses; {b['compounds_three_comparable_doses']} have ≥3. {b['compounds_with_auc']} compounds support eligible AUC, {b['compounds_with_cmax']} support Cmax and {b['compounds_with_both_auc_cmax']} support both somewhere in the validated cohort; {b['both_endpoint_series']} *individual series* support both. AUC/Cmax eligibility is endpoint-specific.

{table(sources,['source','collected_records','eligible_records','compounds','auc_compounds','cmax_compounds','source_publications'])}

Per-source compound counts overlap and must not be summed. Identity uses standardized full InChIKey, with source structures retained; stereochemistry is preserved. Connectivity groups provide a conservative sensitivity count and future split grouping. No names alone are merged. Mixtures and invalid structures are excluded from A and remain in the ledger. Reversible protonation is normalized; original structures and source keys remain available. No tautomer enumeration or largest-fragment guessing is used.

## Verified source findings

- **ChEMBL 37:** public CC BY-SA 3.0 numerical activity API plus original assay/document metadata. Fixed first 10,000 activity IDs from a rat-target ADME query for CL, Vd, T1/2, AUC and Cmax; the query contains {manifest['query_total_records']:,} records. This bounded tranche is neither exhaustive nor random. API filters alone are insufficient: assay organism, source validity, endpoint descriptions, units, censoring and explicit in vivo evidence are independently checked. [API documentation](https://chembl.gitbook.io/chembl-interface-documentation/web-services) and [license](https://chembl.gitbook.io/chembl-interface-documentation/about).
- **PKSmart:** downloaded the publisher's CC BY 4.0 workbook and imported Table S2 only (371 source rows). This table reports raw-scale CL and Vss despite the paper's description of later log-transformed modeling. Units are interpreted from table headers; values are not exponentiated. Human tables, fu-only data and Table S6 predicted animal values are excluded. Missing original dose/publication metadata prevents cross-dose use. [Paper](https://doi.org/10.1186/s13321-025-01066-5), [licensed supplement](https://doi.org/10.6084/m9.figshare.30224986.v1).
- **2026 9,098-compound PBPK study:** internal Roche gRED data shared with pRED; dose-normalized AUC and median repeated endpoints. The listed supplement provides ECCS thresholds, not 9,098 structures/measurements. No compound-level dataset was obtained or counted. [Original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12958287/).
- **NCATS FRDB:** official ZIP successfully downloaded and inspected: 13,455 PK rows, all Homo sapiens. None count toward preclinical A. The archive readme does not establish a blanket redistribution license; no human rows are republished. [Official download](https://drugs.ncats.io/downloads-public).
- **REFRAMEdb:** the public catalog advertises PK annotations but is a clinical-stage library, and catalog PK counts do not prove preclinical numerical availability. The accessible client page did not establish an openly licensed bulk numerical PK export. No annotation count or screening dose-response is imported as PK. [Portal](https://reframedb.org/).
- **CEBS/PMC and PK-DB:** existing reviewed CEBS/PMC endpoints are incorporated once via B; no new CEBS series claimed. PK-DB's documented account/download and reuse restrictions remain unchanged; no retry or bypass was attempted. [Source inventory](data_inventory.md).

## Harmonization and exclusions

{table([dict(reason=k,records=v) for k,v in sorted(reasons.items())],['reason','records'])}

Flags overlap; count excluded rows from the summary rather than summing reasons. Original values/units, standard values/units, descriptions, assay parameters, API IDs, URLs, hashes, retrieval dates and source licenses remain available. Unsupported dimensions, nonpositive/nonfinite values, interval/censored measurements, invalid structures, conflicting identifiers and absent/contradictory in vivo context remain excluded. Source-suspected duplicates and incomplete experimental fields remain flagged. No endpoint is computed from a concentration curve and no missing endpoint is imputed.

Dose-normalization has three explicit states: yes, no (reviewed B), and not_reported. Bodyweight-normalized CL/Vd is not labeled dose-normalized. A dose-normalized exposure never supplies a measured administered dose; no multiplication by a guessed dose is performed. Cmax nM conversions use the source structure molecular weight and are flagged. AUC definition, matrix and route remain distinct. Oral clearance retains an apparent-clearance warning; unspecified Vd/half-life phase is not relabeled Vss/terminal half-life.

A 26-record source-metadata consistency review is preserved in [molecular_description_audit.csv](molecular_description_audit.csv), with expected eligibility decisions in the curation configuration. It identified mixed-route descriptions and effect-time AUC ambiguity; the resulting exclusions apply throughout the tranche. This is not an original-paper validation of all records. The audit metadata retains ChEMBL CC BY-SA 3.0 terms.

## Original-dose search and strict eligibility

The metadata screen found **{len(candidates)} unverified endpoint-specific multi-dose candidates**, written to `dose_candidates_unverified.csv`. Matching publication/compound/species/sex/route/matrix is only a discovery lead: original protocol, formulation, regimen, sampling and AUC definitions still require review. These candidates contribute **zero new validated B series**. Exact-DOI Europe PMC queries for the seven candidate publications found no PMC full text and reported no open-access status; see [access audit](molecular_candidate_access.json). This is not proof that no lawful copy exists elsewhere. Original protocols were not obtained, including conjugate/formulation and AUC-window ambiguities. PKSmart's original 2013 compilation is cited by the paper, but per-record original study/dose mapping is absent from Table S2; it cannot establish dose series.

## Compound-disjoint modeling feasibility

There are {feasibility['scaffold_groups']} scaffold groups (including one shared empty scaffold for {feasibility['acyclic_compounds']} acyclic compounds). The known-publication subset has {feasibility['known_publication_dependency_components']} connected compound–publication components; the largest contains {feasibility['largest_component_compounds']} connectivity groups. These counts describe potential dependency structure, not independent animal experiments.

Large-scale molecular modeling is now structurally possible for sufficiently populated homogeneous endpoint/species/route cohorts. `prediction_cohort_feasibility.csv` separates definitions, units and dose-normalization status. Before freezing evaluation, resolve within-source and cross-source duplicates, review original study metadata and restrict exposure modeling to compatible definitions and known doses. Group all records of a molecular identity, connectivity variants and linked publications together. PKSmart rows with unknown original studies cannot support a study-disjoint evaluation claim. AUC/Cmax at unknown doses cannot be treated as unconditional molecular properties. No models or train/test splits were created in this expansion.

## Reproduce and continue

Run `PYTHONPATH=src python -m pktk.molecular` once online after the existing strict pipeline; use `--offline` for cached replay. Run `PYTHONPATH=src python scripts/verify_molecular_reproducibility.py` to compare A artifacts and verify B numerical artifacts are unchanged. Raw caches and generated data remain ignored by Git. Source-specific licenses apply to derived datasets; do not distribute the combined data under the code MIT license.

Next: adjudicate the dose-candidate queue against original papers; map PKSmart original records; review larger, versioned ChEMBL tranches across additional species. Source study validation takes precedence over increasing the structure count. Dataset B analysis remains in [scientific_analysis.md](scientific_analysis.md).
'''
    (root/'reports/molecular_pk_expansion.md').write_text(report)
    (root/'reports/molecular_qc.md').write_text('# Dataset A QC\n\n'+table([dict(reason=k,records=v) for k,v in sorted(reasons.items())],['reason','records'])+f'\n\nCollected {len(x)} endpoint records; eligible {len(good)}; excluded {len(x)-len(good)}. Flags overlap. Every source endpoint remains in endpoint_records.csv. No B promotion occurs.\n')
    return summary
