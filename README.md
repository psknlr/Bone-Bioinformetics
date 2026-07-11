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

## Detailed analysis workflow

### 0. Inputs and reproducibility controls

| Input | Source | Used for | Output |
|---|---|---|---|
| `osteoporosis_extraction_output_finalV6.xlsx` | Curated workbook in this repository | Classical record filtering, strict de-duplication, herb parsing, historical stability | `deduplicated_classical_records.csv`, `stable_herb_modules.csv`, `historical_stability.csv` |
| Public REST/API responses | PubChem, ChEMBL, Open Targets, STRING, GWAS Catalog | Modern pharmacology, disease genetics and network evidence | cached under `.cache/`, summarized as CSV tables |
| Public expression matrices | PanglaoDB marker table, GSE224152, GSE246769 | Reference marker overlap, marrow expression localisation, osteoclast dynamics | `reference_marker_overlap_localisation.csv`, `gse224152_target_expression.csv`, `gse246769_osteoclast_dynamics.csv` |
| Large public-resource manifest | GEO, LINCS/scPerturb, PrimeKG/BindingDB, PRIDE/ProteomeXchange, MetaboLights | Full-scale planning and optional Colab downloads | `full_scale_resource_execution_plan.csv`, `full_scale_download_log.csv` |

Reproducibility rules:

1. PNG figures are generated at runtime and ignored by git.
2. Large raw files are not committed; full mode records a download/execution plan.
3. High-order analyses requiring unavailable full summary statistics/QTL/perturbation inputs are marked as pending rather than simulated.
4. API responses and downloaded small matrices are cached locally under `.cache/` for repeatability.

### 1. Stable classical herb-module discovery

**Dataset.** `osteoporosis_extraction_output_finalV6.xlsx`, sheet `Sheet1`.

**Processing.**

1. Keep records with `status == ok` and `consider_include == True`.
2. Construct a strict `record_key` from `Books | year | Diagnosis | Chapter | Title` and drop duplicates.
3. Parse `therapeutic_drugs_json`, normalize herb/formula names and build one herb set per record.
4. Assign dynasty strata from the `year` column.

**Methods.**

- FP-Growth itemset mining with minimum support bounded by record count.
- Lift against independence.
- Random permutation null model preserving marginal herb frequencies.
- Benjamini-Hochberg FDR correction.
- Stratified persistence summaries by dynasty, TCM syndrome and diagnosis.

**Primary outputs.**

- `stable_herb_modules.csv`
- `historical_stability.csv`
- `Fig1_stable_modules.png`
- `Fig2_historical_stability.png`

### 2. Real component-target evidence and network medicine

**Core herbs.** `杜仲`, `牛膝`, `续断`, `骨碎补`.

**Compound query seeds.** Literature-guided compound names are used only as search terms; chemical IDs and target evidence are fetched from public databases.

| Herb | Query terms |
|---|---|
| 杜仲 | pinoresinol diglucoside, aucubin, geniposidic acid, chlorogenic acid |
| 牛膝 | 20-hydroxyecdysone, ecdysterone, oleanolic acid, ginsenoside Ro |
| 续断 | asperosaponin VI, akebia saponin D, loganin, chlorogenic acid |
| 骨碎补 | naringin, neoeriocitrin, eriodictyol, kaempferol |

**Data sources.**

| Source | Access mode | Role |
|---|---|---|
| PubChem PUG REST | compound name property endpoint | CID, title, formula, canonical SMILES |
| ChEMBL Web Services | molecule search, activity, target endpoints | ChEMBL molecule IDs, human activity records, pChEMBL, target gene symbols |
| Open Targets Platform GraphQL | disease `MONDO_0005298` | osteoporosis-associated genes and scores |
| STRING API | protein network endpoint | target-disease PPI network proximity |

**Methods.**

1. Query PubChem for compound identity.
2. Query ChEMBL for molecule IDs and activity records.
3. Retain human activity evidence with `pChEMBL >= 5`.
4. Resolve ChEMBL targets to gene symbols via target component synonyms.
5. Query Open Targets for osteoporosis-associated genes.
6. Query STRING for a protein interaction network covering drug targets and osteoporosis genes.
7. Compute shortest-path proximity from herb target sets to the osteoporosis target module.

**Primary outputs.**

- `pubchem_chembl_compounds.csv`
- `component_target_evidence_tiers.csv`
- `opentargets_osteoporosis_targets.csv`
- `network_proximity.csv`
- `Fig3_network_proximity.png`

### 3. Reference marker and real expression localisation

Current localisation is explicitly reported as **reference marker-based cell-type overlap localization** unless it comes from a downloaded expression matrix.

#### 3.1 PanglaoDB marker overlap

**Source.** `https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz`

**Cells retained.** Osteoblasts, osteoclasts, osteoclast precursor cells, osteocytes and stromal cells with human marker support.

**Method.** Intersect ChEMBL-derived target gene symbols with downloaded PanglaoDB markers and the curated bone marker panel.

**Output.** `reference_marker_overlap_localisation.csv`

#### 3.2 GSE224152 marrow non-haematopoietic expression

**Source.** `GSE224152_matrix-Combined-Loic-GenesInRows.csv.gz` from GEO.

**Biological use.** MSC/stromal/vascular localisation in human marrow non-haematopoietic cells.

**Method.** For every candidate target gene present in the matrix, compute mean expression, detection fraction and maximum expression across cells.

**Output.** `gse224152_target_expression.csv`

#### 3.3 GSE246769 osteoclast differentiation dynamics

**Source.** `GSE246769_RNAcounts_Differentiation.txt.gz` from GEO.

**Biological use.** Multi-donor human CD14+ monocyte-to-osteoclast differentiation dynamics across d0, d2, d5 and d9.

**Method.** Parse gene symbols from the annotation column, compute log2(CPM+1) per sample and average target expression by differentiation day; report d9-d0 change.

**Output.** `gse246769_osteoclast_dynamics.csv`

### 4. Genetics-anchored causal pharmacology prioritisation

**Sources.**

| Source | Role |
|---|---|
| GWAS Catalog REST API | bone trait study discovery for osteoporosis, BMD and fracture |
| Open Targets Platform | osteoporosis target association scores |
| ChEMBL | experimental compound-target evidence |
| PanglaoDB/GSE224152/GSE246769 | cell and expression context |

**Methods.**

1. Query GWAS Catalog for osteoporosis, bone mineral density and fracture-related studies.
2. Combine ChEMBL pChEMBL evidence, Open Targets osteoporosis scores, marker overlap, GSE224152 detection and GSE246769 osteoclast dynamics.
3. Assign evidence-count summaries and an A/B/C candidate label:
   - `A_candidate_requires_coloc_MR_confirmation`: strong multi-evidence support plus Open Targets evidence, still requiring fine-mapping/coloc/MR confirmation.
   - `B_multi_omics_partial_support`: genetics or several omics layers partially support the target.
   - `C_pharmacology_only_or_weak_omics`: mainly pharmacology evidence or weak omics support.
4. Estimate preliminary pharmacology direction from osteoclast dynamics/expression heuristics and explicitly mark fine-mapping/coloc/MR as pending full summary-statistics/QTL inputs.

**Outputs.**

- `gwas_catalog_bone_trait_studies.csv`
- `genetics_anchored_causal_pharmacology_scores.csv`
- `human_genetics_prioritised_targets.csv`
- `Fig5_human_genetics.png`
- `Fig6_causal_evidence_heatmap.png`

### 5. Full-scale single-cell, spatial, perturbation and knowledge-graph planning

Run:

```bash
python scripts/nature_bone_pipeline.py \
  --workbook osteoporosis_extraction_output_finalV6.xlsx \
  --outdir results \
  --analysis-scope full
```

This writes a full execution plan without downloading multi-GB files by default. To download optional large resources in Colab, run:

```bash
python scripts/nature_bone_pipeline.py \
  --workbook osteoporosis_extraction_output_finalV6.xlsx \
  --outdir results \
  --analysis-scope full \
  --download-large \
  --large-data-dir data/full_scale
```

**Full-scale resources covered.**

| Module | Public resources |
|---|---|
| Bone marrow and bone single-cell/spatial | GSE253355, GSE224152, GSE147287, GSE147390, GSE169396, GSE255646, GSE242414 |
| Cross-species trajectories | GSE269583, GSE145477 |
| Osteoclast dynamics | GSE246769 |
| Causal genetics | GWAS Catalog, Open Targets credible sets/L2G/colocalisation, eQTL Catalogue |
| Perturbation biology | LINCS L1000 Phase 1/2, scPerturb |
| Knowledge graph and binding | PrimeKG, BindingDB, ChEMBL, PubChem |
| Proteomics/metabolomics validation | PXD035745, PXD017804, MTBLS11650 |

**Outputs.**

- `public_expression_dataset_manifest.csv`
- `high_order_public_resource_manifest.csv`
- `full_scale_resource_execution_plan.csv`
- `full_scale_download_log.csv`
- `Fig7_high_order_resource_map.png`

### 6. Figure generation

The pipeline writes publication-style figures to `results/figures/` at runtime:

| Figure | Content |
|---|---|
| `Fig1_stable_modules.png` | permutation/FDR-supported herb modules |
| `Fig2_historical_stability.png` | historical persistence across dynasties |
| `Fig3_network_proximity.png` | target proximity to osteoporosis network |
| `Fig4_single_cell_localisation.png` | marker-overlap localisation summary |
| `Fig5_human_genetics.png` | Open Targets-supported candidate genes |
| `Fig6_causal_evidence_heatmap.png` | multi-evidence causal pharmacology heatmap |
| `Fig7_high_order_resource_map.png` | public high-order resource network |
| `Fig8_osteoclast_dynamic_targets.png` | GSE246769 osteoclast target dynamics |

PNG files are ignored by git, but the Colab notebook regenerates and displays them automatically.

### 7. Colab reproduction

Click the badge at the top of this README or open:

```text
https://colab.research.google.com/github/pariskang/Bone-Bioinformetics/blob/main/notebooks/bone_bioinformatics_colab.ipynb
```

The notebook:

1. clones `pariskang/Bone-Bioinformetics` if needed;
2. installs `requirements.txt`;
3. runs the full planning pipeline;
4. previews the core CSV tables;
5. displays all generated figures.
