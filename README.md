# Bone-Bioinformetics

Reproducible analysis workflow for integrating classical Chinese medicine records with modern osteoporosis network biology.

## Main workflow

Run the full analysis from the curated workbook:

```bash
python scripts/nature_bone_pipeline.py --workbook osteoporosis_extraction_output_finalV6.xlsx --outdir results
```

The pipeline implements the requested main narrative:

1. **Stable herb modules**: filters included classical records, creates strict record hashes, removes duplicate book/year/diagnosis/chapter/title entries, mines herb modules with FP-Growth, estimates lift against independence, calculates permutation p values, and applies Benjamini-Hochberg FDR correction.
2. **Historical stability**: compares the Du-Zhong/Niu-Xi/Xu-Duan/Gu-Sui-Bu core across dynasty, TCM syndrome, and diagnosis strata.
3. **Component-target and network medicine layer**: exports a target evidence-tier table for the core herbs and compares the core module's osteoporosis-network proximity with single herbs and random target sets.
4. **Cell and human genetics layer**: exports single-cell localisation scores and a prioritised target table designed for downstream GEO/CELLxGENE pseudo-bulk, UCell, trajectory, GWAS enrichment, colocalisation, and cis-MR replacement with project-specific external data.

## Outputs

Generated tables are written to `results/tables/`:

- `deduplicated_classical_records.csv`
- `stable_herb_modules.csv`
- `historical_stability.csv`
- `component_target_evidence_tiers.csv`
- `network_proximity.csv`
- `single_cell_localisation.csv`
- `human_genetics_prioritised_targets.csv`

Generated figures are written to `results/figures/` when the pipeline is run locally. PNG binaries are intentionally ignored and not committed:

- `Fig1_stable_modules.png`
- `Fig2_historical_stability.png`
- `Fig3_network_proximity.png`
- `Fig4_single_cell_localisation.png`
- `Fig5_human_genetics.png`

## Notes for manuscript-grade extension

The workbook is sufficient for reproducible historical module mining. The script also creates transparent evidence templates for HERB/HIT 2.0/ChEMBL/BindingDB/PubChem, GEO/CELLxGENE, and GWAS analyses so those external source tables can be swapped into the same workflow without changing the overall study design.
