# Bone-Bioinformetics

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pariskang/Bone-Bioinformetics/blob/main/notebooks/bone_bioinformatics_colab.ipynb)

Reproducible analysis workflow for integrating classical Chinese medicine records with modern osteoporosis network biology.

## Main workflow

Run the full analysis from the curated workbook:

```bash
python scripts/nature_bone_pipeline.py --workbook osteoporosis_extraction_output_finalV6.xlsx --outdir results
# full manifest/planning mode for all large public datasets
python scripts/nature_bone_pipeline.py --workbook osteoporosis_extraction_output_finalV6.xlsx --outdir results --analysis-scope full
```

The pipeline implements the requested main narrative:

1. **Stable herb modules**: filters included classical records, creates strict record hashes, removes duplicate book/year/diagnosis/chapter/title entries, mines herb modules with FP-Growth, estimates lift against independence, calculates permutation p values, and applies Benjamini-Hochberg FDR correction.
2. **Historical stability**: compares the Du-Zhong/Niu-Xi/Xu-Duan/Gu-Sui-Bu core across dynasty, TCM syndrome, and diagnosis strata.
3. **Component-target and network medicine layer**: queries PubChem/ChEMBL for core-herb compound evidence and compares the core module's osteoporosis-network proximity with single herbs using STRING protein interactions and Open Targets osteoporosis genes.
4. **Cell and human genetics layer**: exports reference marker-based cell-type overlap localisation, downloads the real GSE224152 marrow non-haematopoietic expression matrix for target-expression summaries, downloads the real GSE246769 multi-donor osteoclast differentiation bulk RNA-seq matrix for dynamic validation, and keeps UCell/pseudo-bulk/trajectory reserved for true single-cell matrices with adequate metadata.

## Outputs

Generated tables are written to `results/tables/`:

- `deduplicated_classical_records.csv`
- `stable_herb_modules.csv`
- `historical_stability.csv`
- `pubchem_chembl_compounds.csv`
- `component_target_evidence_tiers.csv`
- `opentargets_osteoporosis_targets.csv`
- `network_proximity.csv`
- `single_cell_localisation.csv`
- `reference_marker_overlap_localisation.csv`
- `gse224152_target_expression.csv`
- `gse246769_osteoclast_dynamics.csv`
- `public_expression_dataset_manifest.csv`
- `human_genetics_prioritised_targets.csv`
- `gwas_catalog_bone_trait_studies.csv`
- `genetics_anchored_causal_pharmacology_scores.csv`
- `high_order_public_resource_manifest.csv`
- `full_scale_resource_execution_plan.csv`
- `full_scale_download_log.csv`

Generated figures are written to `results/figures/` when the pipeline is run locally. PNG binaries are intentionally ignored and not committed:

- `Fig1_stable_modules.png`
- `Fig2_historical_stability.png`
- `Fig3_network_proximity.png`
- `Fig4_single_cell_localisation.png`
- `Fig5_human_genetics.png`
- `Fig6_causal_evidence_heatmap.png`
- `Fig7_high_order_resource_map.png`
- `Fig8_osteoclast_dynamic_targets.png`

## Notes for manuscript-grade extension

The workbook is sufficient for reproducible historical module mining. Modern pharmacology and human-genetics layers are now fetched from public APIs instead of simulated data: PubChem and ChEMBL for compound/target evidence, Open Targets for osteoporosis-associated human genes, and STRING for protein-network proximity. If an API returns no evidence for a compound, the output records the absence rather than fabricating targets.


## Expression-data policy

Current cell-localisation outputs are named **reference marker-based cell-type overlap localization** unless they come from a downloaded expression matrix. The workflow downloads real GSE224152 and GSE246769 matrices for expression summaries, records a manifest of recommended GEO datasets, and explicitly marks GSE249471 as bulk-only/not suitable for single-cell UCell, pseudo-bulk, or trajectory analyses.

## High-order modules

The workflow now adds a first-pass high-order evidence layer using public data only: GWAS Catalog trait-study discovery for BMD/osteoporosis/fracture, a genetics-anchored causal pharmacology scoring table that combines ChEMBL activity, Open Targets osteoporosis evidence, marker overlap, GSE224152 expression and GSE246769 osteoclast dynamics, and a manifest for larger follow-up resources such as eQTL Catalogue, LINCS L1000, scPerturb, GSE284089, PrimeKG, BindingDB and proteomics datasets. Fine-mapping, coloc, SMR/HEIDI and MR are explicitly marked as pending full summary statistics/QTL inputs rather than simulated.

## Reproducibility and audit

A reproducible Colab notebook is provided at `notebooks/bone_bioinformatics_colab.ipynb`. The notebook installs dependencies, runs the real-data pipeline, displays generated figures, and previews the key CSV outputs. See `docs/REAL_DATA_AUDIT.md` for a source-by-source audit of which analyses use real data and which high-order modules are intentionally marked as pending full summary-statistics/QTL/perturbation inputs rather than simulated.

## Full-scale mode

Use `--analysis-scope full` to emit a complete execution plan for all recommended large public resources, including large single-cell/spatial matrices, LINCS/scPerturb, PrimeKG/BindingDB and proteomics/metabolomics validation datasets. By default this mode records the plan without downloading multi-GB resources. Add `--download-large --large-data-dir data/full_scale` in Colab only when you intentionally want to fetch the optional large files.

The plotting layer now generates additional Nature-style bioinformatics figures at runtime: a multi-evidence causal score heatmap, a high-order public-resource graph and osteoclast differentiation target dynamics. PNG outputs remain ignored by git and are regenerated in Colab.
