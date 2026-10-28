"""Generate current reports from actual processed data, retaining scope caveats."""
import json
import pandas as pd


def expanded_report(root,s,series,attrition,resolved):
    primary=series[series.primary_eligible]
    table=['| Compound | Species / sex | Route | Doses | Eligible endpoints | Study |',
           '|---|---|---|---:|---|---|']
    for _,r in primary.iterrows():
        ep=', '.join(e for e,ok in [('AUC',r.eligible_auc),('Cmax',r.eligible_cmax)] if ok)
        table.append(f'| {r.compound_name} | {r.species} / {r.sex} | {r.route} | {r.n_doses} | {ep} | {r.study_id} |')
    stages=json.loads((root/'reports/stage_snapshots.json').read_text())
    stage_table=['| Scope after stage | Validated compounds | Independent studies | Series | AUC | Cmax | Both eligible |',
                 '|---|---:|---:|---:|---:|---:|---:|']
    for v in stages:
        stage_table.append('| '+' | '.join(str(v[k]) for k in ['stage','validated_compounds','independent_studies','validated_series','auc_series','cmax_series','both_eligible_series'])+' |')
    text=f'''# Verified expanded dataset — 2026-10-28

The primary subset now contains **{s['primary_comparable_series']} comparable multi-dose series, {s['primary_comparable_compound_identifiers']} compound identities and {s['primary_independent_studies']} independent source studies**. The previous run had 9 series, 7 compounds and 5 studies. Species, sex, route and repeated-dose day strata within a study are correlated; they are not independent studies.

## Counts after each stage

{chr(10).join(stage_table)}

These are reproducible scope snapshots using the final identity registry, not model results. Source studies are identified by publication DOI for PMC and experiment number for CEBS; a CEBS source and its CvTdb copy are one study.

- Collected parent identifiers: **{s['preclinical_parent_compound_identifiers']}**; resolved unique full-InChIKey structures: **{s['resolved_unique_parent_structures']}**; unresolved identifiers: **{s['unresolved_parent_identifiers']}**. A fully resolved unique-molecule total cannot be claimed. Analyte identifiers including metabolites: {s['preclinical_analyte_identifiers_collected']}.
- Collected document-level study proxies: **{s['preclinical_document_studies']}**; this is not the number of independent experiments. PMC publications collected: {s['collected_pmc_publications']}.
- Compounds with ≥2 comparable doses: **{s['compounds_at_least_two_comparable_doses']}**; with ≥3: **{s['compounds_at_least_three_comparable_doses']}**.
- AUC-eligible series: **{s['primary_auc_eligible_series']}**; Cmax-eligible: **{s['primary_cmax_eligible_series']}**; both eligible: **{s['primary_both_endpoints_eligible_series']}**. {s['primary_series_with_auc_and_cmax']} primary series report both endpoints, including AUC definitions not accepted for comparison.
- Primary dose groups: **{s['primary_endpoint_groups']}/{s['curated_endpoint_groups']} ({100*s['fraction_curated_groups_primary']:.1f}%)**. Verified comparable parent-identifier fraction: **{100*s['fraction_parent_identifiers_verified_comparable']:.2f}%**. These denominators are dose groups and parent identifiers, not raw observations.
- Collected species: {', '.join(s['collected_species'])}; routes: {', '.join(s['collected_routes'])}. Primary subset: {', '.join(s['primary_species'])}; {', '.join(s['primary_routes'])}; plasma.

## Validated series

{chr(10).join(table)}

## Source contributions and low-yield approaches

CEBS/NTP group spreadsheets add 14 validated series from three already collected compounds and three source studies: benzophenone (4 AUC series), formamide (8 AUC series) and isoeugenol (2 mouse Cmax series). These validate existing CvTdb candidates; they are not new independent copies of those observations. Formamide's predicted Cmax values remain in original source cells and are not presented as observed Cmax.

Two additional PMC publications add three series: lorlatinib (rat, oral; AUC0-24h and Cmax) and methyl eugenol (mouse, oral; day 1 AUC0-inf/Cmax and day 28 Cmax). The two methyl eugenol time-course strata use the same animals and count as one study. Day 28 AUC0-t is preserved but not eligible without an explicit endpoint window. Lorlatinib is CC BY-NC 4.0; its local derived data retain the noncommercial restriction and are not committed. It can be excluded using the license flag.

The composite audit reviewed 18/112 candidates across three protocols: 14 supported group-endpoint comparisons, two had plasma/blood conflicts and two had 35-versus-70 mg/kg dose conflicts. The remaining 94 were not reviewed. This purposive sample is not an estimate of yield for all 112 candidates. No unreviewed group is promoted and no individual-animal trajectory is invented. See [audit decisions](composite_audit.csv).

The identity pass does not itself create experiments. Direct PubChem REST remains outside the allowed robots paths; PubChem CIDs are exact-structure UniChem cross-references. EPA DSSTox bulk matches require exact DTXSID/CAS and full source-SMILES/InChIKey agreement. Cached ChEBI evidence is checked for cross-source conflicts; remaining identifiers use ChEBI queries. Unsuccessful and conflicting candidates remain explicit. PK-DB previously yielded metadata but no eligible preclinical measurements; the DA-6886 supplement added no exposure series. The simmitecan XML endpoint failed after retries in this expansion and contributed no data. Additional search leads are not counted as collected studies.

## QC and remaining limitations

The previously unqueried {s.get('previously_unqueried_identifiers','pending stage recount')} identifiers now have statuses: `{json.dumps(s.get('previously_unqueried_status_counts',{}),sort_keys=True)}`. Full registry identity status counts: `{json.dumps(s['identity_status_counts'],sort_keys=True)}`. {s['resolved_pubchem_cids']} registry entries have one PubChem CID and {s['primary_compounds_with_resolved_structures']} primary compounds have resolved structures. Full stereochemical/isotopic keys are retained; salts and ambiguous mappings are not collapsed. CID lookups cover primary/candidate compounds and existing caches; other resolved compounds may have no optional CID query. See [EPA bulk manifest](dsstox_identity_manifest.json).

The pipeline retains all {s['cvtdb_original_observations']:,} original CvTdb observations, including {s['preclinical_observations']:,} preclinical observations. The original audit found {attrition['singleton_series']:,} singleton series out of {attrition['source_series_ids']:,}; sparse sampling and censored observations explain much of the curve attrition. Endpoint imports preserve all original spreadsheet fields/cells. CEBS actual versus nominal dose, repeated model fits, predicted versus observed peaks, and source metadata errors are explicitly distinguished.

{ s['flagged_protocol_series'] } curated series remain outside the primary subset. Primary records still flag missing fasting information, unavailable AUC extrapolated fractions, composite/aggregate uncertainty and source exceptions. Different AUC definitions, study days and matrices are never pooled. Missing structures and identifier conflicts remain limitations. The 112-candidate universe was exploratory, not 112 established experimental series.

No complex models were trained. The cohort has only {s['primary_comparable_compound_identifiers']} compounds in {s['primary_independent_studies']} studies; several compounds share a publication. This is insufficient for meaningful compound-disjoint tuning or scaffold evaluation. Existing baseline results and splits concern only the original HY-071085 pilot. Expanded endpoints have not been used to claim new model performance.

See [inventory](data_inventory.md), [QC](qc_report.md), [validation](validation.md), [statistics](expansion_summary.json), [identity provenance](identity_source_manifest.json) and [unresolved issues](unresolved_issues.json).
'''
    (root/'reports/final_report.md').write_text(text)
    (root/'reports/collection_summary.md').write_text('# Current collection summary\n\nSee [final report](final_report.md) for stage counts, source contributions and limitations.\n\n```json\n'+json.dumps(s,indent=2)+'\n```\n')
    qc_path=root/'reports/qc_counts.json';qc=json.loads(qc_path.read_text())
    qc['expansion']={k:s[k] for k in ['curated_endpoint_groups','primary_endpoint_groups','primary_comparable_series','flagged_protocol_series','extractions_failed','source_overlap_groups','identity_status_counts','reviewed_composite_candidates','validated_composite_series']}
    qc_path.write_text(json.dumps(qc,indent=2)+'\n')
    (root/'reports/qc_report.md').write_text('# Current QC report\n\nNo source records are silently deleted. Eligibility applies to curated endpoints, not removal from source caches.\n\n```json\n'+json.dumps(qc,indent=2)+'\n```\n\nExact source/table checksums and protocol phrases protect reviewed mappings. All 112 exploratory groups have an audit status; 18 were reviewed. Group endpoints are not individual measurements. Dose/matrix conflicts are withheld. Formamide Cmax_pred is not Cmax_obs. Benzophenone uses actual dose comments, not nominal dose. Isoeugenol rat dose conflicts prevent primary eligibility. AUC units and definitions are endpoint-specific; 0-t without a verified terminal time is not compared. Source exceptions, license terms, original cells and predicted values are retained.\n')
    (root/'reports/benchmark_feasibility.md').write_text('# Current benchmark feasibility\n\n'+chr(10).join(stage_table)+f'\n\nOnly {s["primary_comparable_compound_identifiers"]} compounds and {s["primary_independent_studies"]} independent source studies are available. No complex ML training or compound/scaffold split is justified. Task A is structurally feasible for the eligible endpoint-specific groups, but expanded episodes have not been evaluated. Existing splits and baselines are HY-071085 pilot-only. Task D lacks independently validated nonlinear labels; group means and source variability cannot establish mechanism.\n\nFuture evaluation must group full molecular identity and related source publications together, keep repeat-dose days and sex/route strata correlated, separate observed Cmax from predicted Cmax, and respect noncommercial source terms.\n')
    issue_path=root/'reports/unresolved_issues.json'
    issues=json.loads(issue_path.read_text())
    updates={
        'protocols':('partially_resolved','18 of 112 candidates audited against three CEBS source protocols: 14 validated using group endpoints, two matrix conflicts, two dose conflicts; 94 unreviewed.'),
        'structures':('partially_resolved',f"Full registry assessed: {s['resolved_unique_parent_structures']} resolved unique full-key structures, {s['unresolved_parent_identifiers']} unresolved identifiers. Exact EPA DTXSID/CAS and source-key validation; cached ChEBI evidence compared. Nonstandard source InChI options are not guessed."),
        'auc':('partially_resolved',f"{s['primary_auc_eligible_series']} primary series have comparable AUC definitions. Unspecified terminal windows and unknown tail fractions remain limitations. Predicted Cmax is excluded from observed-Cmax endpoints."),
        'sample_size':('open',f"{s['primary_comparable_compound_identifiers']} primary compounds across {s['primary_independent_studies']} source studies remain insufficient for meaningful compound-disjoint molecular ML. Related strata/publications require grouping."),
        'other_sources':('partially_resolved','CEBS adds 14 validated series and two new PMC publications add three series. Simmitecan XML failed after retries; rotundic acid XML collected but CC BY-NC-ND terms require separate review. FDA/EMA wider extraction remains future work.'),
        'overlap':('partially_resolved',f"{s['source_overlap_groups']} selected source overlap groups identified by document/DOI and identity. CEBS endpoints validate existing CvTdb studies and are not independent replicates. Remaining incomplete identifiers/protocol links prevent global deduplication proof."),
        'source_quality':('open','Isoeugenol rat middle dose is 35 mg/kg in one field and 70 in group labels/CvTdb. Benzophenone plasma/blood conflicts remain withheld. Historical vehicle/unit ambiguities remain. Composite endpoints preserve source methods and uncertainty.'),
        'distribution_terms':('documented','EPA DSSTox CC0; CRAN package GPL-3; source PMC CC BY or CC BY-NC. Lorlatinib derivatives remain noncommercial, raw/generated data are not committed. Rotundic acid not imported pending CC BY-NC-ND review.'),
        'identity_access':('documented','Direct PubChem REST and ChEBI FTP paths remain disallowed by robots. Public EPA DSSTox bulk CSV, ChEBI and UniChem APIs are used. CIDs are optional and not queried for every resolved identifier.')}
    for issue in issues:
        if issue['id'] in updates:issue['status'],issue['detail']=updates[issue['id']]
    issue_path.write_text(json.dumps(issues,indent=2)+'\n')
