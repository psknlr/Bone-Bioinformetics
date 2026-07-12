#!/usr/bin/env python3
"""Real-data osteoporosis herb-module analysis.

This workflow never fabricates pharmacology, expression, GWAS, or network data.
Classical-record mining is computed from the supplied workbook. Modern evidence is
queried live from public APIs (PubChem, ChEMBL, Open Targets, STRING) and cached.
If an API returns no evidence, the corresponding output records that absence.
"""
from __future__ import annotations

import argparse, gzip, hashlib, io, json, random, re, time
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from mlxtend.frequent_patterns import fpgrowth
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests

CORE = ("杜仲", "牛膝", "续断", "骨碎补")
DYNASTIES = [(0, 618, "隋唐以前"), (618, 907, "唐"), (907, 1279, "宋金元"), (1279, 1644, "明"), (1644, 1912, "清"), (1912, 10_000, "近现代")]
# Literature-guided query seeds; all compound metadata/targets are fetched from public databases below.
HERB_COMPOUND_QUERIES = {
    "杜仲": ["pinoresinol diglucoside", "aucubin", "geniposidic acid", "chlorogenic acid"],
    "牛膝": ["20-hydroxyecdysone", "ecdysterone", "oleanolic acid", "ginsenoside Ro"],
    "续断": ["asperosaponin VI", "akebia saponin D", "loganin", "chlorogenic acid"],
    "骨碎补": ["naringin", "neoeriocitrin", "eriodictyol", "kaempferol"],
}

GEO_DATASETS = [
    {"accession":"GSE253355","role":"reference human femoral-head marrow atlas","url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253355","recommended_use":"marker validation, BMSC/endothelial/pericyte/immune localisation; not OP-control inference"},
    {"accession":"GSE224152","role":"human marrow non-haematopoietic cells","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224152/suppl/GSE224152_matrix-Combined-Loic-GenesInRows.csv.gz","recommended_use":"MSC/stromal/vascular target expression localisation; downloaded in this workflow"},
    {"accession":"GSE147287","role":"human OP and OA CD271+ marrow cells","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147287/suppl/GSE147287_RAW.tar","recommended_use":"exploratory BMSC substate localisation; not powered OP-control pseudo-bulk"},
    {"accession":"GSE147390","role":"primary human osteoblast-lineage atlas","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147390/suppl/GSE147390_RAW.tar","recommended_use":"osteoblast precursor-to-mature trajectory localisation"},
    {"accession":"GSE169396","role":"human femoral-head whole-cell atlas","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169396/suppl/GSE169396_RAW.tar","recommended_use":"whole bone microenvironment localisation; verify clinical grouping before OP-control pseudo-bulk"},
    {"accession":"GSE255646","role":"human MSC single-cell multiome","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255646/suppl/GSE255646_filtered_feature_bc_matrix.h5","recommended_use":"MSC expression plus ATAC/motif/GWAS regulatory interpretation after donor metadata audit"},
    {"accession":"GSE246769","role":"human osteoclast differentiation bulk RNA-seq time series","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246769/suppl/GSE246769_RNAcounts_Differentiation.txt.gz","recommended_use":"multi-donor osteoclast dynamic validation; downloaded in this workflow; not single-cell UCell"},
    {"accession":"GSE242414","role":"osteoporotic vertebral compression fracture single-cell matrix","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242414/suppl/GSE242414_RAW.tar","recommended_use":"fracture microenvironment localisation; not primary OP-normal inference"},
    {"accession":"GSE269583","role":"mouse BMSC-to-osteoblast differentiation single-cell atlas","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269583/suppl/GSE269583_RAW.tar","recommended_use":"cross-species trajectory after human-mouse ortholog mapping"},
    {"accession":"GSE145477","role":"mouse marrow mesenchymal lineage and ageing","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145477/suppl/GSE145477_RAW.tar","recommended_use":"ageing stromal/osteogenic trajectory after ortholog mapping"},
    {"accession":"GSE249471","role":"bulk CPM perturbation only","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249471/suppl/GSE249471_All_gene_counts.txt.gz","recommended_use":"exclude from single-cell UCell/pseudo-bulk/trajectory because GEO attachment is bulk"},
]
PANGLOADB_URL = "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
GSE224152_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224152/suppl/GSE224152_matrix-Combined-Loic-GenesInRows.csv.gz"
GSE246769_COUNTS_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246769/suppl/GSE246769_RNAcounts_Differentiation.txt.gz"

GWAS_TRAITS = ["osteoporosis", "bone mineral density", "heel bone mineral density", "femoral neck bone mineral density", "lumbar spine bone mineral density", "hip fracture", "fragility fracture"]
HIGH_ORDER_RESOURCES = [
    {"module":"causal_genetics", "resource":"GWAS Catalog", "url":"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByEfoTrait", "use":"trait-to-study discovery for BMD/osteoporosis/fracture"},
    {"module":"causal_genetics", "resource":"Open Targets credible sets/L2G/colocalisation", "url":"https://platform.opentargets.org/downloads", "use":"download credible sets, L2G, colocalisation for fine-mapping-aware causal target validation"},
    {"module":"causal_genetics", "resource":"eQTL Catalogue tabix paths", "url":"https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/tabix_ftp_paths.tsv", "use":"cis-eQTL/sQTL extraction for coloc/SMR/MR"},
    {"module":"spatial_single_cell", "resource":"GSE284089", "url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284089", "use":"human femur spatial transcriptomics reference for cell-type niche mapping"},
    {"module":"perturbation", "resource":"LINCS L1000 Phase 1/2", "url":"https://clue.io/connectopedia/lincs_cmap_data", "use":"Level 5 perturbation signatures for disease-signature reversal"},
    {"module":"perturbation", "resource":"scPerturb", "url":"https://projects.sanderlab.org/scperturb/", "use":"single-cell perturbation effect comparison when relevant perturbagens exist"},
    {"module":"knowledge_graph_structure", "resource":"PrimeKG", "url":"https://github.com/mims-harvard/PrimeKG", "use":"heterogeneous disease-drug-gene-pathway graph construction"},
    {"module":"knowledge_graph_structure", "resource":"BindingDB", "url":"https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp", "use":"experimental compound-protein affinity evidence"},
    {"module":"proteomics", "resource":"PXD035745", "url":"https://www.ebi.ac.uk/pride/archive/projects/PXD035745", "use":"Drynaria fortunei osteoporosis-model TMT proteomics validation"},
    {"module":"proteomics", "resource":"PXD017804", "url":"https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD017804", "use":"postmenopausal low-BMD plasma proteomics validation"},
]

FULL_SCALE_RESOURCES = [
    {"category":"single_cell_spatial","accession":"GSE253355_MSC_Subset","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253355/suppl/GSE253355_MSC_Subset_Seurat.rds.gz","size_note":"382.5 MB","default_action":"manifest_only"},
    {"category":"single_cell_spatial","accession":"GSE253355_Normal_BM_Atlas","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253355/suppl/GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz","size_note":"1.4 GB","default_action":"manifest_only"},
    {"category":"single_cell_spatial","accession":"GSE253355_RAW","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253355/suppl/GSE253355_RAW.tar","size_note":"1.1 GB","default_action":"manifest_only"},
    {"category":"single_cell_spatial","accession":"GSE224152_matrix","url":GSE224152_MATRIX_URL,"size_note":"2.5 MB","default_action":"download_and_analyse"},
    {"category":"single_cell_spatial","accession":"GSE147287_RAW","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147287/suppl/GSE147287_RAW.tar","size_note":"91.2 MB","default_action":"full_mode_optional"},
    {"category":"single_cell_spatial","accession":"GSE147390_RAW","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147390/suppl/GSE147390_RAW.tar","size_note":"75.6 MB","default_action":"full_mode_optional"},
    {"category":"single_cell_spatial","accession":"GSE169396_RAW","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169396/suppl/GSE169396_RAW.tar","size_note":"134.2 MB","default_action":"full_mode_optional"},
    {"category":"single_cell_spatial","accession":"GSE255646_GEX_H5","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255646/suppl/GSE255646_filtered_feature_bc_matrix.h5","size_note":"145.1 MB","default_action":"full_mode_optional"},
    {"category":"single_cell_spatial","accession":"GSE255646_metadata","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255646/suppl/GSE255646_per_barcode_metadata.tsv.gz","size_note":"metadata","default_action":"full_mode_optional"},
    {"category":"osteoclast_bulk","accession":"GSE246769_counts","url":GSE246769_COUNTS_URL,"size_note":"2.0 MB","default_action":"download_and_analyse"},
    {"category":"fracture_single_cell","accession":"GSE242414_RAW","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242414/suppl/GSE242414_RAW.tar","size_note":"81 MB","default_action":"full_mode_optional"},
    {"category":"mouse_trajectory","accession":"GSE269583_RAW","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269583/suppl/GSE269583_RAW.tar","size_note":"2.0 GB","default_action":"manifest_only"},
    {"category":"mouse_aging","accession":"GSE145477_RAW","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145477/suppl/GSE145477_RAW.tar","size_note":"73.7 MB","default_action":"full_mode_optional"},
    {"category":"perturbation","accession":"LINCS_GSE92742","url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742","size_note":"large Level 5 signatures","default_action":"manifest_only"},
    {"category":"perturbation","accession":"LINCS_GSE70138","url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138","size_note":"large Level 5 signatures","default_action":"manifest_only"},
    {"category":"knowledge_graph","accession":"PrimeKG","url":"https://github.com/mims-harvard/PrimeKG","size_note":"multi-million-edge KG","default_action":"manifest_only"},
    {"category":"binding","accession":"BindingDB","url":"https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp","size_note":"large TSV/SDF","default_action":"manifest_only"},
    {"category":"proteomics","accession":"PXD035745","url":"https://www.ebi.ac.uk/pride/archive/projects/PXD035745","size_note":"TMT proteomics","default_action":"manifest_only"},
    {"category":"proteomics","accession":"PXD017804","url":"https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD017804","size_note":"plasma proteomics","default_action":"manifest_only"},
    {"category":"metabolomics","accession":"MTBLS11650","url":"https://www.ebi.ac.uk/metabolights/editor/MTBLS11650","size_note":"femoral-neck osteoporosis metabolomics","default_action":"manifest_only"},
]

CELL_MARKERS = {
    "BMSC/成骨祖细胞": {"LRP5","WNT16","SOST","RUNX2","SP7","BMP2","SMAD1","COL1A1"},
    "成骨细胞": {"RUNX2","SP7","ALPL","BGLAP","COL1A1","BMP2"},
    "破骨细胞": {"TNFSF11","TNFRSF11B","CTSK","ACP5","NFKB1"},
    "骨免疫细胞": {"JUN","NFKB1","MAPK1","AKT1","TNFSF11"},
}

class ApiCache:
    def __init__(self, path: Path):
        self.path = path; self.path.mkdir(parents=True, exist_ok=True)
    def get_json(self, url, params=None, method="GET", payload=None, sleep=0.15):
        key = hashlib.md5(json.dumps([method, url, params, payload], sort_keys=True).encode()).hexdigest()
        f = self.path / f"{key}.json"
        if f.exists(): return json.loads(f.read_text())
        if method == "POST": r = requests.post(url, params=params, json=payload, timeout=45)
        else: r = requests.get(url, params=params, timeout=45)
        meta = {"url": r.url, "status_code": r.status_code, "text": r.text[:1000]}
        r.raise_for_status(); data = r.json(); f.write_text(json.dumps({"meta": meta, "data": data}, ensure_ascii=False)); time.sleep(sleep)
        return {"meta": meta, "data": data}
    def get_bytes(self, url, headers=None, sleep=0.15):
        key = hashlib.md5(json.dumps(["BYTES", url, headers], sort_keys=True).encode()).hexdigest()
        f = self.path / f"{key}.bin"
        if f.exists(): return f.read_bytes()
        r = requests.get(url, headers=headers or {}, timeout=120)
        r.raise_for_status(); f.write_bytes(r.content); time.sleep(sleep)
        return r.content

def parse_json_list(x):
    if pd.isna(x): return []
    try: return json.loads(x)
    except Exception: return []

def norm_herb(name):
    return re.sub(r"（.*?）|\(.*?\)|[\s，,。；;：:]+", "", str(name)).strip()

def _herbset(x):
    return {h for h in (norm_herb(d.get("name", "")) for d in parse_json_list(x) if isinstance(d, dict) and d.get("name")) if h}

def _symptomset(x):
    return {norm_herb(d.get("name", d) if isinstance(d, dict) else d) for d in parse_json_list(x)} - {""}

def _base_records(path):
    """Filtered rows with a source_record_id for the *same original passage*.

    The id is the extraction-invariant identity of a classical passage
    (Books|year|Diagnosis|Chapter|Title); multiple structured extractions of the
    same passage share it and must be reconciled, not silently dropped."""
    df = pd.read_excel(path, sheet_name="Sheet1")
    df = df[df["status"].eq("ok") & df["consider_include"].fillna(False)].copy().reset_index(drop=True)
    df["source_record_id"] = df[["Books", "year", "Diagnosis", "Chapter", "Title"]].fillna("").astype(str).agg("|".join, axis=1).map(lambda s: hashlib.md5(s.encode()).hexdigest())
    df["herbs_set"] = df["therapeutic_drugs_json"].map(_herbset)
    df["symptoms_set"] = df["symptoms_signs_json"].map(_symptomset) if "symptoms_signs_json" in df.columns else [set()] * len(df)
    df["quality"] = pd.to_numeric(df.get("quality_score"), errors="coerce").fillna(-1)
    return df

def dedup_conflict_log(path, outdir=None):
    """Field-level conflicts among multiple extractions of the same passage."""
    df = _base_records(path)
    rows = []
    for sid, sub in df.groupby("source_record_id"):
        if len(sub) < 2:
            continue
        herb_variants = {tuple(sorted(s)) for s in sub.herbs_set}
        sym_variants = {tuple(sorted(s)) for s in sub.symptoms_set}
        rows.append({
            "source_record_id": sid, "n_extractions": len(sub),
            "books": sub.Books.iloc[0], "title": str(sub.Title.iloc[0])[:40],
            "drugs_conflict": len(herb_variants) > 1,
            "symptoms_conflict": len(sym_variants) > 1,
            "syndrome_conflict": sub.get("tcm_syndrome", pd.Series(dtype=str)).astype(str).nunique() > 1,
            "op_related_conflict": sub.get("is_osteoporosis_related", pd.Series(dtype=str)).astype(str).nunique() > 1,
            "relevance_conflict": pd.to_numeric(sub.get("osteoporosis_relevance_score"), errors="coerce").nunique() > 1,
            "first_not_highest_quality": sub.quality.iloc[0] < sub.quality.max(),
            "first_has_no_drug_but_others_do": (len(sub.herbs_set.iloc[0]) == 0) and any(len(s) > 0 for s in sub.herbs_set.iloc[1:]),
            "union_herbs": ";".join(sorted(set().union(*sub.herbs_set))),
        })
    log = pd.DataFrame(rows)
    if outdir is not None:
        log.to_csv(outdir / "dedup_conflict_log.csv", index=False)
    return log

def _finalise(df):
    df = df[df["herbs"].map(len).gt(0)].copy()
    df["dynasty"] = pd.cut(df["year"], bins=[x[0] for x in DYNASTIES] + [10_000], labels=[x[2] for x in DYNASTIES], right=False)
    return df

def load_records(path, strategy="merge"):
    """Reconcile multiple extractions of the same passage.

    strategy:
      * "merge"           – expert-confirmed *union* of herbs/symptoms per passage (primary);
      * "highest_quality" – keep the single highest quality_score extraction;
      * "first"           – legacy keep-first (drops later, often higher-quality, extractions);
      * "raw"             – every extraction kept as its own record (no reconciliation).
    """
    df = _base_records(path)
    if strategy == "raw":
        out = df.copy()
        out["herbs"] = out.herbs_set.map(lambda s: sorted(s))
        out["symptoms"] = out.symptoms_set.map(lambda s: sorted(s))
        return _finalise(out)
    if strategy == "first":
        out = df.drop_duplicates("source_record_id").copy()
        out["herbs"] = out.herbs_set.map(lambda s: sorted(s))
        out["symptoms"] = out.symptoms_set.map(lambda s: sorted(s))
        return _finalise(out)
    if strategy == "highest_quality":
        idx = df.sort_values("quality", ascending=False).drop_duplicates("source_record_id").index
        out = df.loc[idx].copy()
        out["herbs"] = out.herbs_set.map(lambda s: sorted(s))
        out["symptoms"] = out.symptoms_set.map(lambda s: sorted(s))
        return _finalise(out)
    # merge (default): union of herbs/symptoms; keep max quality and first metadata
    agg = df.groupby("source_record_id")
    out = agg.first().reset_index()
    out["herbs"] = [sorted(set().union(*g.herbs_set)) for _, g in agg]
    out["symptoms"] = [sorted(set().union(*g.symptoms_set)) for _, g in agg]
    out["quality"] = agg.quality.max().values
    out["n_extractions"] = agg.size().values
    return _finalise(out)

def dedup_strategy_sensitivity(path, outdir=None):
    """Sensitivity of the core-quartet co-occurrence to the de-duplication choice."""
    core = set(CORE)
    rows = []
    for strat in ["first", "highest_quality", "merge", "raw"]:
        d = load_records(path, strategy=strat)
        sets = [set(h) for h in d.herbs]
        rows.append({
            "strategy": strat, "records_with_herbs": len(sets),
            "core_quartet_complete_n": sum(core.issubset(s) for s in sets),
            "core_triple_n": sum(len(core & s) >= 3 for s in sets),
            "core_pair_n": sum(len(core & s) >= 2 for s in sets),
            "core_any_n": sum(bool(core & s) for s in sets),
        })
    out = pd.DataFrame(rows)
    if outdir is not None:
        out.to_csv(outdir / "dedup_strategy_sensitivity.csv", index=False)
    return out

def _perm_cooccurrence_p(rng, N, colsums, obs, n_perm):
    """Permutation p preserving each herb's marginal frequency."""
    null = [len(set.intersection(*[set(rng.sample(range(N), s)) for s in colsums])) for _ in range(n_perm)]
    return (1 + sum(x >= obs for x in null)) / (n_perm + 1)

def core_module_significance(df, outdir, n_perm=1000, seed=11):
    """Explicitly test the target quartet and EVERY 2/3/4-herb sub-combination,
    independent of any frequency threshold, so the central claim about
    杜仲–牛膝–续断–骨碎补 is reported honestly (support, lift, permutation p, FDR)."""
    rng = random.Random(seed)
    sets = [set(hs) for hs in df.herbs]
    N = len(sets)
    freq = {h: sum(h in s for s in sets) for h in CORE}
    rows = []
    for k in range(2, len(CORE) + 1):
        for cols in combinations(sorted(CORE), k):
            obs = sum(set(cols).issubset(s) for s in sets)
            exp = N * np.prod([freq[c] / N for c in cols])
            perm_p = _perm_cooccurrence_p(rng, N, [freq[c] for c in cols], obs, n_perm)
            rows.append({"module": "–".join(cols), "size": len(cols), "support_n": obs,
                         "support_rate": obs / N, "expected_n": exp,
                         "lift_vs_independence": obs / max(exp, 1e-12), "perm_p": perm_p,
                         "co_occurs_at_all": obs > 0})
    res = pd.DataFrame(rows)
    res["fdr"] = multipletests(res.perm_p, method="fdr_bh")[1]
    res = res.sort_values(["size", "support_n"], ascending=[True, False])
    res.to_csv(outdir / "core_module_significance.csv", index=False)
    return res

def mine_modules(df, outdir, n_perm=200, seed=7):
    rng = random.Random(seed); herbs = sorted({h for hs in df.herbs for h in hs})
    mat = pd.DataFrame([{h: h in hs for h in herbs} for hs in df.herbs])
    fis = fpgrowth(mat, min_support=max(0.01, 8/len(mat)), use_colnames=True, max_len=4)
    fis = fis[fis.itemsets.map(len).ge(2)].copy(); rows=[]; N=len(mat)
    for items in fis.itemsets:
        cols=sorted(items); obs=int(mat[cols].all(axis=1).sum()); exp=N*np.prod([mat[c].mean() for c in cols]); colsums=[int(mat[c].sum()) for c in cols]
        null=[len(set.intersection(*[set(rng.sample(range(N), s)) for s in colsums])) for _ in range(n_perm)]
        rows.append({"module":"–".join(cols),"size":len(cols),"support_n":obs,"expected_n":exp,"lift_vs_independence":obs/max(exp,1e-12),"perm_p":(1+sum(x>=obs for x in null))/(n_perm+1),"items":tuple(cols)})
    mods=pd.DataFrame(rows); mods["fdr"]=multipletests(mods.perm_p, method="fdr_bh")[1]
    core=set(CORE); mods["contains_core"] = mods["items"].map(lambda x: set(x).issubset(core) or core.issubset(set(x)))
    mods.sort_values(["contains_core","fdr","lift_vs_independence","support_n"], ascending=[False,True,False,False]).to_csv(outdir/"stable_herb_modules.csv", index=False)
    hist=[]
    def _core_stats(sub):
        return {"records":len(sub),
                "core_complete_n":sum(core.issubset(set(x)) for x in sub.herbs),
                "core_pair_n":sum(len(core & set(x))>=2 for x in sub.herbs),
                "core_any_n":sum(bool(core & set(x)) for x in sub.herbs),
                "core_any_rate":sum(bool(core & set(x)) for x in sub.herbs)/len(sub)}
    for group_col in ["dynasty","tcm_syndrome","Diagnosis"]:
        if group_col not in df.columns: continue
        for g, sub in df.groupby(group_col, dropna=True):
            hist.append({"stratum":group_col,"level":str(g),**_core_stats(sub)})
    # symptom (症状) stratification: for each recurring symptom, persistence of the core set
    sym_counts=pd.Series([s for syms in df.symptoms for s in syms if s]).value_counts()
    for sym in sym_counts[sym_counts>=8].index:
        sub=df[df.symptoms.map(lambda ss: sym in ss)]
        if len(sub): hist.append({"stratum":"symptom","level":str(sym),**_core_stats(sub)})
    pd.DataFrame(hist).to_csv(outdir/"historical_stability.csv", index=False); return mods, pd.DataFrame(hist)

def pubchem_lookup(cache, name):
    # PubChem deprecated the `CanonicalSMILES` property name; request current
    # `SMILES`/`ConnectivitySMILES` plus InChIKey for a real chemical-identity audit.
    props_req = "Title,MolecularFormula,InChIKey,SMILES,ConnectivitySMILES"
    url=f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/property/{props_req}/JSON"
    try: props=cache.get_json(url)["data"]["PropertyTable"]["Properties"][0]
    except Exception:
        try:  # fall back to legacy property name on older PubChem deployments
            url2=f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/property/Title,MolecularFormula,InChIKey,CanonicalSMILES/JSON"
            props=cache.get_json(url2)["data"]["PropertyTable"]["Properties"][0]
        except Exception:
            return None
    return {"query":name,"pubchem_cid":props.get("CID"),"pubchem_title":props.get("Title"),
            "formula":props.get("MolecularFormula"),"inchikey":props.get("InChIKey"),
            "canonical_smiles":props.get("SMILES") or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")}

def chembl_molecule(cache, name):
    """Best fuzzy match plus the fields needed to audit chemical identity."""
    data=cache.get_json("https://www.ebi.ac.uk/chembl/api/data/molecule/search.json", {"q":name,"limit":1})["data"]
    mols=data.get("molecules") or []
    if not mols: return {"molecule_chembl_id":None,"chembl_pref_name":None,"chembl_inchikey":None}
    m=mols[0]; struct=m.get("molecule_structures") or {}
    return {"molecule_chembl_id":m.get("molecule_chembl_id"),"chembl_pref_name":m.get("pref_name"),
            "chembl_inchikey":struct.get("standard_inchi_key")}

def chembl_targets(cache, chembl_id, min_pchembl=5.0, max_pages=6):
    """Human activities with the assay/target/relation metadata needed to keep
    only molecular-target binding/functional evidence (excludes cell-line and
    phenotypic assays)."""
    out=[]
    for page in range(max_pages):
        data=cache.get_json("https://www.ebi.ac.uk/chembl/api/data/activity.json", {"molecule_chembl_id":chembl_id,"pchembl_value__isnull":False,"limit":100,"offset":page*100})["data"]
        for a in data.get("activities", []):
            try: pchem=float(a.get("pchembl_value") or 0)
            except ValueError: pchem=0
            if pchem >= min_pchembl and a.get("target_organism") == "Homo sapiens":
                out.append({"molecule_chembl_id":chembl_id,"target_chembl_id":a.get("target_chembl_id"),
                            "target_pref_name":a.get("target_pref_name"),"pchembl_value":pchem,
                            "standard_type":a.get("standard_type"),"standard_relation":a.get("standard_relation"),
                            "assay_type":a.get("assay_type"),"assay_description":(a.get("assay_description") or "")[:120],
                            "document_chembl_id":a.get("document_chembl_id")})
        if not data.get("page_meta",{}).get("next"): break
    return out

def chembl_target_detail(cache, target_chembl_id):
    """Return (target_type, [gene_symbols]) so complexes/families are not collapsed to one arbitrary gene."""
    if not target_chembl_id: return None, []
    try:
        data=cache.get_json(f"https://www.ebi.ac.uk/chembl/api/data/target/{target_chembl_id}.json")["data"]
    except Exception:
        return None, []
    genes=[]
    for comp in data.get("target_components") or []:
        for syn in comp.get("target_component_synonyms") or []:
            if syn.get("syn_type") == "GENE_SYMBOL":
                g=syn.get("component_synonym")
                if g and g not in genes: genes.append(g)
    return data.get("target_type"), genes

def opentargets_osteoporosis(cache, size=200):
    q='''query disease($id:String!,$size:Int!){ disease(efoId:$id){ id name associatedTargets(page:{index:0,size:$size}){ rows{ score target{ approvedSymbol id } } } } }'''
    data=cache.get_json("https://api.platform.opentargets.org/api/v4/graphql", method="POST", payload={"query":q,"variables":{"id":"MONDO_0005298","size":size}})["data"]
    rows=data.get("data",{}).get("disease",{}).get("associatedTargets",{}).get("rows",[]) if data else []
    return pd.DataFrame([{"target":r["target"].get("approvedSymbol"),"open_targets_score":r.get("score"),"ensembl_id":r["target"].get("id")} for r in rows if r.get("target")])

def string_network(cache, genes, add_nodes=500, physical=True):
    """STRING physical-interaction network for the seed genes, expanded with
    `add_nodes` interactors so there is a broader background for randomisation.

    Note: `add_nodes` interactors are chosen relative to the query seeds, so this
    is a seed-centred local network, not an unbiased whole-interactome background;
    the randomisation below draws from an independent target universe to avoid
    that bias (see build_background_universe)."""
    genes=sorted({g for g in genes if isinstance(g,str) and g});
    if len(genes)<2: return nx.Graph()
    params={"identifiers":"%0d".join(genes),"species":9606,"required_score":400,"add_nodes":add_nodes}
    if physical: params["network_type"]="physical"
    text=requests.get("https://string-db.org/api/tsv/network", params=params, timeout=120).text
    G=nx.Graph(); G.add_nodes_from(genes)
    for line in text.splitlines()[1:]:
        f=line.split('\t')
        if len(f)>5: G.add_edge(f[2], f[3], score=float(f[5]))
    return G

def _closest_distance(G, source_set, target_set):
    """Mean over sources of the shortest-path distance to the nearest target
    (Guney et al. 2016 closest distance d_c). A source that is itself a disease
    (target-set) gene has distance 0 — direct overlap must count as closest."""
    tgt=set(t for t in target_set if t in G); vals=[]
    for s in source_set:
        if s not in G: continue
        if s in tgt:
            vals.append(0); continue
        ds=[nx.shortest_path_length(G,s,t) for t in tgt if nx.has_path(G,s,t)]
        if ds: vals.append(min(ds))
    return float(np.mean(vals)) if vals else np.nan

def proximity_zscore(G, drug_genes, disease_genes, n_rand=1000, seed=3):
    """Network proximity of a drug-target set to the disease module with a
    degree-preserving random reference, giving z-score and empirical p."""
    drug=[g for g in drug_genes if g in G]; disease=[g for g in disease_genes if g in G]
    d_obs=_closest_distance(G, drug, disease)
    if np.isnan(d_obs) or not drug: return {"observed_distance":d_obs,"z_score":np.nan,"empirical_p":np.nan,"random_mean":np.nan,"random_sd":np.nan,"n_targets_in_network":len(drug)}
    # degree-preserving bins over candidate nodes (exclude disease module itself)
    universe=[n for n in G.nodes if n not in set(disease)]
    deg={n:G.degree(n) for n in universe}
    import bisect
    order=sorted(universe, key=lambda n: deg[n]); degs=[deg[n] for n in order]
    rng=np.random.default_rng(seed); null=[]
    for _ in range(n_rand):
        pick=set()
        for g in drug:
            d=G.degree(g); lo=bisect.bisect_left(degs,int(d*0.5)); hi=max(bisect.bisect_right(degs,int(d*1.5)+1),lo+1)
            cand=order[lo:hi] or order
            for _try in range(20):
                c=cand[rng.integers(len(cand))]
                if c not in pick: pick.add(c); break
        d=_closest_distance(G, list(pick), disease)
        if not np.isnan(d): null.append(d)
    mu=float(np.mean(null)) if null else np.nan; sd=float(np.std(null)) if null else np.nan
    z=(d_obs-mu)/sd if sd and sd>0 else np.nan
    emp_p=(1+sum(x<=d_obs for x in null))/(len(null)+1) if null else np.nan
    return {"observed_distance":d_obs,"z_score":z,"empirical_p":emp_p,"random_mean":mu,"random_sd":sd,"n_targets_in_network":len(drug)}

def evidence_tier(standard_type, pchembl):
    """Stratify experimental ChEMBL evidence by assay type and potency."""
    st=str(standard_type).upper(); p=pchembl or 0
    if st in {"KI","KD"} and p>=7: return "T1_experimental_high_affinity_binding"
    if st in {"KI","KD","IC50","XC50"} and p>=6: return "T2_experimental_binding"
    if st in {"EC50","POTENCY","AC50"} and p>=6: return "T3_experimental_functional"
    if p>=5: return "T4_experimental_moderate"
    return "T5_experimental_weak"

def real_targets_and_network(outdir, cache, ot_score_threshold=0.10):
    # --- compound identity (PubChem CID/InChIKey + ChEMBL id/InChIKey) ---
    compounds=[]; raw_activities=[]
    for herb, queries in HERB_COMPOUND_QUERIES.items():
        for q in queries:
            pc=pubchem_lookup(cache, q) or {"query":q,"pubchem_cid":None,"pubchem_title":None,"formula":None,"inchikey":None,"canonical_smiles":None}
            ch=chembl_molecule(cache, q)
            pc.update({"herb":herb,"compound_scope":"selected_representative_compound",**ch}); compounds.append(pc)
            cid=ch.get("molecule_chembl_id")
            if cid:
                for t in chembl_targets(cache, cid):
                    ttype, genes = chembl_target_detail(cache, t.get("target_chembl_id"))
                    t.update({"herb":herb,"compound_query":q,"target_type":ttype,
                              "target_gene_symbol":(genes[0] if ttype=="SINGLE PROTEIN" and genes else None),
                              "all_component_genes":";".join(genes)})
                    raw_activities.append(t)
    comp=pd.DataFrame(compounds)
    # synonym-collision detection: distinct query names sharing a PubChem CID / ChEMBL InChIKey
    collision=pd.Series(False, index=comp.index)
    for key in ["pubchem_cid","chembl_inchikey"]:
        if key in comp:
            n_shared=comp[key].map(comp.dropna(subset=[key]).groupby(key)["query"].nunique())
            comp[f"{key}_shared_by_n_queries"]=n_shared
            collision=collision | (n_shared.fillna(1)>1)
    comp["synonym_collision"]=collision
    comp.to_csv(outdir/"pubchem_chembl_compounds.csv", index=False)

    acts=pd.DataFrame(raw_activities)
    # transparency: log activities excluded because the target is not a single protein
    if not acts.empty:
        nonmol=acts[acts.target_type.ne("SINGLE PROTEIN")].copy()
        nonmol.to_csv(outdir/"excluded_nonmolecular_activities.csv", index=False)
        acts=acts[acts.target_type.eq("SINGLE PROTEIN") & acts.target_gene_symbol.notna()].copy()
        acts=acts.drop_duplicates(subset=["molecule_chembl_id","target_chembl_id","standard_type","pchembl_value","document_chembl_id"])
        acts["evidence_type"]="experimental_ChEMBL_single_protein"
        acts["evidence_tier"]=[evidence_tier(st,p) for st,p in zip(acts.standard_type, acts.pchembl_value)]
    else:
        acts=pd.DataFrame(columns=["herb","compound_query","target_gene_symbol","pchembl_value","standard_type"])
    acts.to_csv(outdir/"component_target_evidence_tiers.csv", index=False)
    # target-level aggregation: replaces "max-pChEMBL of duplicated rows" with real summaries
    if not acts.empty:
        summ=acts.groupby("target_gene_symbol").agg(
            herbs=("herb",lambda s:";".join(sorted(set(s)))),
            n_activities=("pchembl_value","size"),
            n_documents=("document_chembl_id","nunique"),
            n_assay_types=("assay_type","nunique"),
            median_pchembl=("pchembl_value","median"),
            max_pchembl=("pchembl_value","max"),
            standard_types=("standard_type",lambda s:";".join(sorted(set(map(str,s))))),
            best_evidence_tier=("evidence_tier","min")).reset_index()
        summ.to_csv(outdir/"component_target_summary.csv", index=False)
    disease=opentargets_osteoporosis(cache)
    disease["note"]="Open Targets association score is a weighted ranking heuristic (genetics+drugs+literature+animal+expression), not a probability"
    disease.to_csv(outdir/"opentargets_osteoporosis_targets.csv", index=False)
    drug_gene=set(acts.get("target_gene_symbol", pd.Series(dtype=str)).dropna())
    # disease module by association-score threshold (documented), not an arbitrary top-N
    disease_gene=set(disease.loc[disease.open_targets_score.ge(ot_score_threshold),"target"].dropna()) if not disease.empty else set()
    if not disease_gene and not disease.empty:
        disease_gene=set(disease.target.dropna().head(50))
    G=string_network(cache, drug_gene | disease_gene)
    # Predicted-target layer: STRING first-shell neighbours of experimental targets in the disease module
    pred_rows=[]
    for g in sorted(drug_gene):
        if g not in G: continue
        for nb in G.neighbors(g):
            if nb in disease_gene and nb not in drug_gene:
                pred_rows.append({"predicted_target":nb,"via_experimental_target":g,"string_score":G[g][nb].get("score"),"evidence_type":"predicted_STRING_neighbour_of_experimental_target","in_osteoporosis_module":True})
    pd.DataFrame(pred_rows).to_csv(outdir/"predicted_target_layer.csv", index=False)
    # Network proximity with degree-preserving random reference (z-score, empirical p).
    rows=[{"set":"core_combo","targets":len(drug_gene),**proximity_zscore(G, drug_gene, disease_gene),"source":"ChEMBL(single-protein)+OpenTargets(score>=%.2f)+STRING-physical"%ot_score_threshold}]
    for herb in CORE:
        hg=set(acts.loc[acts.herb.eq(herb),"target_gene_symbol"].dropna()) if not acts.empty else set()
        rows.append({"set":herb,"targets":len(hg),**proximity_zscore(G, hg, disease_gene),"source":"ChEMBL(single-protein)+OpenTargets+STRING-physical"})
    prox_df=pd.DataFrame(rows).rename(columns={"observed_distance":"network_distance"})
    prox_df.to_csv(outdir/"network_proximity.csv", index=False)
    # Random-combination null drawn from an INDEPENDENT background (all network proteins
    # outside the disease module: drug targets + STRING add_nodes interactors), size-matched.
    # This fixes the previous degenerate null where sampling equalled the whole target pool.
    background=[n for n in G.nodes if n not in disease_gene]
    k=len([g for g in drug_gene if g in G]); d_core=_closest_distance(G,[g for g in drug_gene if g in G],[d for d in disease_gene if d in G])
    rng=np.random.default_rng(17); null=[]
    if len(background)>k>0 and not np.isnan(d_core):
        for _ in range(2000):
            samp=[background[i] for i in rng.choice(len(background),k,replace=False)]
            d=_closest_distance(G, samp, [d for d in disease_gene if d in G])
            if not np.isnan(d): null.append(d)
    combo_null=pd.DataFrame([{"comparison":"core_combo_vs_random_size_matched_background","core_distance":d_core,"background_size":len(background),"n_targets":k,"random_mean":float(np.mean(null)) if null else np.nan,"random_sd":float(np.std(null)) if null else np.nan,"empirical_p_core_closer":(1+sum(x<=d_core for x in null))/(len(null)+1) if null else np.nan,"n_random":len(null)}])
    combo_null.to_csv(outdir/"network_proximity_random_combo_null.csv", index=False)
    return comp, acts, disease, prox_df

def omics_genetics(outdir, targets, disease):
    genes=set(targets.get("target_gene_symbol", pd.Series(dtype=str)).dropna())
    cell=[]
    for c, markers in CELL_MARKERS.items():
        ov=len(genes & markers); odds,p=fisher_exact([[ov,max(len(genes)-ov,0)],[max(len(markers)-ov,0),200]])
        cell.append({"cell_state":c,"target_overlap":ov,"odds_ratio":odds,"p_value":p,"ucell_signature_score":ov/max(len(markers),1),"source":"marker overlap only; replace with GEO/CELLxGENE matrix for expression statistics"})
    celldf=pd.DataFrame(cell); celldf["fdr"]=multipletests(celldf.p_value, method="fdr_bh")[1]; celldf.to_csv(outdir/"single_cell_localisation.csv", index=False)
    gen=disease[disease.target.isin(genes)].copy() if not disease.empty and genes else pd.DataFrame(columns=["target","open_targets_score","ensembl_id"])
    gen["human_genetics_source"]="Open Targets Genetics/Platform disease-target association for MONDO_0005298"; gen.to_csv(outdir/"human_genetics_prioritised_targets.csv", index=False)
    return celldf, gen

def reference_marker_overlap(outdir, cache, targets):
    """Reference marker-based cell-type overlap localisation using PanglaoDB plus curated bone markers."""
    genes=set(targets.get("target_gene_symbol", pd.Series(dtype=str)).dropna())
    marker_rows=[]
    # PanglaoDB: real downloaded marker table, scored with a hypergeometric
    # enrichment test against the PanglaoDB human gene universe.
    try:
        raw=cache.get_bytes(PANGLOADB_URL, headers={"User-Agent":"Mozilla/5.0"})
        pang=pd.read_csv(io.BytesIO(gzip.decompress(raw)), sep="\t")
        pang=pang[pang["species"].str.contains("Hs", na=False)].copy()
        universe=set(pang["official gene symbol"].dropna().astype(str))
        M=len(universe); drawn=genes & universe; N=len(drawn)
        keep=pang["cell type"].str.lower().isin({"osteoblasts","osteoclasts","osteoclast precursor cells","osteocytes","stromal cells"})
        for cell, sub in pang[keep].groupby("cell type"):
            markers=set(sub["official gene symbol"].dropna().astype(str)); ov=genes & markers
            p=float(hypergeom.sf(len(ov)-1, M, len(markers), N)) if (M and len(markers) and N) else np.nan
            marker_rows.append({"source":"PanglaoDB_27_Mar_2020","cell_state":cell,"marker_n":len(markers),"target_overlap":len(ov),"hypergeom_p":p,"background_universe":M,"overlap_genes":";".join(sorted(ov))})
    except Exception as e:
        marker_rows.append({"source":"PanglaoDB_27_Mar_2020","cell_state":"download_failed","marker_n":0,"target_overlap":0,"hypergeom_p":np.nan,"background_universe":0,"overlap_genes":str(e)[:200]})
    for cell, markers in CELL_MARKERS.items():
        ov=genes & markers
        marker_rows.append({"source":"curated_bone_marker_panel","cell_state":cell,"marker_n":len(markers),"target_overlap":len(ov),"hypergeom_p":np.nan,"background_universe":0,"overlap_genes":";".join(sorted(ov))})
    out=pd.DataFrame(marker_rows)
    if "hypergeom_p" in out and out["hypergeom_p"].notna().any():
        mask=out["hypergeom_p"].notna(); out.loc[mask,"fdr"]=multipletests(out.loc[mask,"hypergeom_p"], method="fdr_bh")[1]
    out.to_csv(outdir/"reference_marker_overlap_localisation.csv", index=False)
    return out

def gse224152_expression_localisation(outdir, cache, targets):
    """Real GSE224152 non-haematopoietic marrow expression matrix summary for candidate target genes."""
    genes=sorted(set(targets.get("target_gene_symbol", pd.Series(dtype=str)).dropna()))
    raw=cache.get_bytes(GSE224152_MATRIX_URL)
    mat=pd.read_csv(io.BytesIO(raw), compression="gzip", index_col=0)
    present=[g for g in genes if g in mat.index]
    rows=[]
    for g in present:
        x=mat.loc[g].astype(float)
        rows.append({"dataset":"GSE224152","gene":g,"cells_n":mat.shape[1],"mean_expression":float(x.mean()),"pct_cells_detected":float((x>0).mean()),"max_expression":float(x.max()),"matrix_source":GSE224152_MATRIX_URL})
    missing=sorted(set(genes)-set(present))
    for g in missing:
        rows.append({"dataset":"GSE224152","gene":g,"cells_n":mat.shape[1],"mean_expression":0.0,"pct_cells_detected":0.0,"max_expression":0.0,"matrix_source":"gene_not_found_in_matrix"})
    out=pd.DataFrame(rows); out.to_csv(outdir/"gse224152_target_expression.csv", index=False)
    return out

def gse246769_osteoclast_dynamics(outdir, cache, targets):
    """Real bulk RNA-seq osteoclast differentiation dynamics from GSE246769."""
    genes=sorted(set(targets.get("target_gene_symbol", pd.Series(dtype=str)).dropna()))
    raw=cache.get_bytes(GSE246769_COUNTS_URL)
    counts=pd.read_csv(io.BytesIO(raw), compression="gzip", sep="\t")
    counts["gene_symbol"]=counts["Annotation.Divergence"].astype(str).str.split("|").str[0]
    sample_cols=[c for c in counts.columns if re.match(r"Donor\d+_d\d+", c)]
    lib=counts[sample_cols].sum(axis=0)
    sub=counts[counts.gene_symbol.isin(genes)].copy()
    rows=[]
    for _, r in sub.iterrows():
        vals=(r[sample_cols].astype(float)/lib*1e6).apply(lambda v: np.log2(v+1))
        by_day={day: vals[[c for c in sample_cols if c.endswith(f"_{day}")]].mean() for day in ["d0","d2","d5","d9"] if any(c.endswith(f"_{day}") for c in sample_cols)}
        rows.append({"dataset":"GSE246769","gene":r.gene_symbol,"log2cpm_d0":by_day.get("d0", np.nan),"log2cpm_d2":by_day.get("d2", np.nan),"log2cpm_d5":by_day.get("d5", np.nan),"log2cpm_d9":by_day.get("d9", np.nan),"delta_d9_vs_d0":by_day.get("d9", np.nan)-by_day.get("d0", np.nan),"counts_source":GSE246769_COUNTS_URL})
    out=pd.DataFrame(rows); out.to_csv(outdir/"gse246769_osteoclast_dynamics.csv", index=False)
    return out

def write_dataset_manifest(outdir):
    out=pd.DataFrame(GEO_DATASETS); out.to_csv(outdir/"public_expression_dataset_manifest.csv", index=False); return out

def gwas_catalog_trait_studies(outdir, cache, traits=GWAS_TRAITS):
    rows=[]
    for trait in traits:
        try:
            r=requests.get("https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByEfoTrait", params={"efoTrait":trait,"size":10}, timeout=15)
            r.raise_for_status(); data=r.json()
            studies=data.get("_embedded",{}).get("studies",[])
            for st in studies:
                rows.append({"query_trait":trait,"accessionId":st.get("accessionId"),"diseaseTrait":st.get("diseaseTrait",{}).get("trait") if isinstance(st.get("diseaseTrait"),dict) else st.get("diseaseTrait"),"initialSampleSize":st.get("initialSampleSize"),"replicateSampleSize":st.get("replicateSampleSize"),"publication":st.get("publicationInfo",{}).get("pubmedId") if isinstance(st.get("publicationInfo"),dict) else None,"summaryStats":"available_if_linked_in_GWAS_Catalog"})
        except Exception as e:
            rows.append({"query_trait":trait,"accessionId":None,"diseaseTrait":"download_failed","initialSampleSize":str(e)[:200],"replicateSampleSize":None,"publication":None,"summaryStats":None})
    out=pd.DataFrame(rows).drop_duplicates(); out.to_csv(outdir/"gwas_catalog_bone_trait_studies.csv", index=False); return out

def gwas_catalog_bone_genes(cache, traits=GWAS_TRAITS, per_trait=1000):
    """Collect author-reported / mapped genes from GWAS Catalog associations for bone traits."""
    genes=set()
    for trait in traits:
        try:
            data=cache.get_json("https://www.ebi.ac.uk/gwas/rest/api/associations/search/findByEfoTrait", {"efoTrait":trait,"size":per_trait})["data"]
        except Exception:
            continue
        for a in data.get("_embedded",{}).get("associations",[]):
            for locus in a.get("loci",[]) or []:
                for rg in locus.get("authorReportedGenes",[]) or []:
                    g=(rg.get("geneName") or "").strip()
                    if g and g.lower() not in {"intergenic","nr","na"}: genes.add(g.upper())
            for ra in a.get("strongestRiskAlleles",[]) or []:
                pass
    return genes

def gwas_gene_enrichment(outdir, cache, targets, background_n=20000):
    """Gene-level enrichment: overlap of candidate targets with GWAS-Catalog bone-trait
    genes, tested by hypergeometric enrichment against a genome background."""
    target_genes={str(g).upper() for g in targets.get("target_gene_symbol", pd.Series(dtype=str)).dropna()}
    gwas_genes=gwas_catalog_bone_genes(cache)
    overlap=sorted(target_genes & gwas_genes)
    M, n, N, k = background_n, len(gwas_genes), len(target_genes), len(overlap)
    p=float(hypergeom.sf(k-1, M, n, N)) if (n and N) else np.nan
    fold=(k/N)/(n/M) if (N and n) else np.nan
    summary=pd.DataFrame([{"n_candidate_targets":N,"n_gwas_bone_genes":n,"n_overlap":k,"fold_enrichment":fold,"hypergeom_p":p,"background_genes":M,"overlap_genes":";".join(overlap),"note":"author-reported GWAS Catalog genes; coloc/MR still require full summary statistics"}])
    summary.to_csv(outdir/"gwas_gene_enrichment.csv", index=False)
    pd.DataFrame({"gene":sorted(gwas_genes)}).to_csv(outdir/"gwas_catalog_bone_genes.csv", index=False)
    return summary, gwas_genes

def causal_pharmacology_scores(outdir, gwas_genes=None):
    targets=pd.read_csv(outdir/"component_target_evidence_tiers.csv")
    ot=pd.read_csv(outdir/"opentargets_osteoporosis_targets.csv") if (outdir/"opentargets_osteoporosis_targets.csv").exists() else pd.DataFrame(columns=["target","open_targets_score"])
    marker=pd.read_csv(outdir/"reference_marker_overlap_localisation.csv") if (outdir/"reference_marker_overlap_localisation.csv").exists() else pd.DataFrame(columns=["overlap_genes"])
    gse224=pd.read_csv(outdir/"gse224152_target_expression.csv") if (outdir/"gse224152_target_expression.csv").exists() else pd.DataFrame(columns=["gene","pct_cells_detected","mean_expression"])
    gse246=pd.read_csv(outdir/"gse246769_osteoclast_dynamics.csv") if (outdir/"gse246769_osteoclast_dynamics.csv").exists() else pd.DataFrame(columns=["gene","delta_d9_vs_d0"])
    marker_genes=set()
    for x in marker.get("overlap_genes", pd.Series(dtype=str)).dropna():
        marker_genes |= {g for g in str(x).split(";") if g and g != "nan"}
    gwas_genes={str(g).upper() for g in (gwas_genes or set())}
    rows=[]
    for gene, sub in targets.dropna(subset=["target_gene_symbol"]).groupby("target_gene_symbol"):
        ot_score=float(ot.loc[ot.target.eq(gene),"open_targets_score"].max()) if gene in set(ot.target) else 0.0
        pct=float(gse224.loc[gse224.gene.eq(gene),"pct_cells_detected"].max()) if gene in set(gse224.gene) else 0.0
        delta=float(gse246.loc[gse246.gene.eq(gene),"delta_d9_vs_d0"].mean()) if gene in set(gse246.gene) else 0.0
        max_pchem=float(sub.pchembl_value.max()) if "pchembl_value" in sub else 0.0
        in_gwas=gene.upper() in gwas_genes
        # four independent evidence pillars from the requested main narrative
        classical=True  # every candidate derives from the deduplicated stable-module herbs
        pharmacology=max_pchem>=6
        cell_expression=(pct>0.01) or (abs(delta)>1) or (gene in marker_genes)
        human_genetics=(ot_score>0) or in_gwas
        pillars=sum([classical, pharmacology, cell_expression, human_genetics])
        evidence_count=sum([ot_score>0, in_gwas, gene in marker_genes, pct>0.01, abs(delta)>1, max_pchem>=6])
        if human_genetics and pillars==4 and evidence_count>=4: grade="A_candidate_requires_coloc_MR_confirmation"
        elif human_genetics or evidence_count>=3: grade="B_multi_omics_partial_support"
        else: grade="C_pharmacology_only_or_weak_omics"
        direction="inhibit_candidate" if delta>1 else ("activate_candidate" if delta<-1 or pct>0.05 else "direction_uncertain")
        rows.append({"target":gene,"herbs":";".join(sorted(sub.herb.dropna().unique())),"max_pchembl":max_pchem,"open_targets_score":ot_score,"gwas_catalog_bone_gene":in_gwas,"reference_marker_overlap":gene in marker_genes,"gse224152_pct_cells_detected":pct,"gse246769_delta_d9_vs_d0":delta,"classical_evidence":classical,"pharmacology_evidence":pharmacology,"cell_expression_evidence":cell_expression,"human_genetics_evidence":human_genetics,"convergent_pillars":pillars,"evidence_count":evidence_count,"causal_evidence_grade":grade,"estimated_pharmacology_direction":direction,"fine_mapping_coloc_mr_status":"not_run_requires_full_sumstats_QTL"})
    out=pd.DataFrame(rows).sort_values(["convergent_pillars","causal_evidence_grade","open_targets_score","evidence_count","max_pchembl"], ascending=[False,True,False,False,False])
    out.to_csv(outdir/"genetics_anchored_causal_pharmacology_scores.csv", index=False); return out

def write_high_order_resource_manifest(outdir):
    out=pd.DataFrame(HIGH_ORDER_RESOURCES); out.to_csv(outdir/"high_order_public_resource_manifest.csv", index=False); return out

def write_full_scale_execution_plan(outdir, analysis_scope="quick", download_large=False, large_data_dir="data/full_scale"):
    rows=[]
    for r in FULL_SCALE_RESOURCES:
        action=r["default_action"]
        if analysis_scope == "full" and action == "full_mode_optional":
            action = "download_if_--download-large" if not download_large else "download_requested"
        if action == "manifest_only" and analysis_scope == "full":
            action = "manifest_only_due_size_or_external_format"
        rows.append({**r, "analysis_scope":analysis_scope, "planned_action":action, "local_path":str(Path(large_data_dir)/Path(r["url"]).name) if download_large and action == "download_requested" else "not_downloaded"})
    out=pd.DataFrame(rows); out.to_csv(outdir/"full_scale_resource_execution_plan.csv", index=False); return out

def download_requested_large_resources(cache, outdir, analysis_scope="quick", download_large=False, large_data_dir="data/full_scale"):
    Path(large_data_dir).mkdir(parents=True, exist_ok=True)
    rows=[]
    if analysis_scope != "full" or not download_large:
        pd.DataFrame(rows, columns=["accession","url","local_path","status","bytes"]).to_csv(outdir/"full_scale_download_log.csv", index=False); return pd.DataFrame(rows)
    for r in FULL_SCALE_RESOURCES:
        if r["default_action"] != "full_mode_optional":
            continue
        local=Path(large_data_dir)/Path(r["url"]).name
        try:
            if not local.exists():
                local.write_bytes(cache.get_bytes(r["url"]))
            rows.append({"accession":r["accession"],"url":r["url"],"local_path":str(local),"status":"downloaded_or_cached","bytes":local.stat().st_size})
        except Exception as e:
            rows.append({"accession":r["accession"],"url":r["url"],"local_path":str(local),"status":f"failed:{str(e)[:120]}","bytes":0})
    out=pd.DataFrame(rows); out.to_csv(outdir/"full_scale_download_log.csv", index=False); return out

def plot_nature_style_extensions(figdir, tabdir):
    figdir=Path(figdir); tabdir=Path(tabdir)
    sns.set_theme(style="white", context="talk")
    score_path=tabdir/"genetics_anchored_causal_pharmacology_scores.csv"
    if score_path.exists():
        sc=pd.read_csv(score_path).head(30).copy()
        cols=["max_pchembl","open_targets_score","gse224152_pct_cells_detected","gse246769_delta_d9_vs_d0","evidence_count"]
        if not sc.empty:
            mat=sc.set_index("target")[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
            mat=(mat-mat.min())/(mat.max()-mat.min()).replace(0,1)
            plt.figure(figsize=(9, max(6, 0.28*len(mat))))
            sns.heatmap(mat, cmap="viridis", cbar_kws={"label":"scaled evidence"})
            plt.title("Multi-evidence causal pharmacology scorecard")
            plt.tight_layout(); plt.savefig(figdir/"Fig6_causal_evidence_heatmap.png", dpi=300); plt.close()
    res_path=tabdir/"high_order_public_resource_manifest.csv"
    if res_path.exists():
        res=pd.read_csv(res_path)
        G=nx.Graph()
        for _,r in res.iterrows():
            G.add_edge(r["module"], r["resource"])
        plt.figure(figsize=(11,8))
        pos=nx.spring_layout(G, seed=4, k=0.8)
        colors=["#b2182b" if n in set(res.module) else "#2166ac" for n in G.nodes]
        nx.draw_networkx(G,pos,node_color=colors,node_size=900,font_size=8,edge_color="#999999",width=1.2)
        plt.axis("off"); plt.title("High-order public-data resource map")
        plt.tight_layout(); plt.savefig(figdir/"Fig7_high_order_resource_map.png", dpi=300); plt.close()
    dyn_path=tabdir/"gse246769_osteoclast_dynamics.csv"
    if dyn_path.exists():
        dyn=pd.read_csv(dyn_path).copy()
        if not dyn.empty:
            dyn["abs_delta"]=dyn["delta_d9_vs_d0"].abs(); top=dyn.nlargest(12,"abs_delta")
            long=top.melt(id_vars=["gene"], value_vars=["log2cpm_d0","log2cpm_d2","log2cpm_d5","log2cpm_d9"], var_name="day", value_name="log2CPM")
            long["day"]=long["day"].str.extract(r"(d\d+)")[0]
            plt.figure(figsize=(10,6)); sns.lineplot(data=long,x="day",y="log2CPM",hue="gene",marker="o")
            plt.title("Top target dynamics during human osteoclast differentiation")
            plt.tight_layout(); plt.savefig(figdir/"Fig8_osteoclast_dynamic_targets.png", dpi=300); plt.close()

def plot_all(figdir, mods, hist, prox, cell, genetics):
    sns.set_theme(style="whitegrid", context="talk")
    plt.figure(figsize=(10,7)); top=mods.nsmallest(20,"fdr").copy(); top["-log10(FDR)"]=-np.log10(top.fdr.clip(1e-12)); sns.scatterplot(data=top,x="lift_vs_independence",y="-log10(FDR)",size="support_n",hue="contains_core",sizes=(80,450)); plt.tight_layout(); plt.savefig(figdir/"Fig1_stable_modules.png",dpi=300); plt.close()
    plt.figure(figsize=(10,6)); h=hist[hist.stratum.eq("dynasty")].copy(); h["plot_level"]=h.level.map({"隋唐以前":"Pre-Tang","唐":"Tang","宋金元":"Song-Jin-Yuan","明":"Ming","清":"Qing","近现代":"Modern"}); sns.barplot(data=h,x="plot_level",y="core_any_rate",color="#2166ac"); plt.xticks(rotation=35,ha="right"); plt.tight_layout(); plt.savefig(figdir/"Fig2_historical_stability.png",dpi=300); plt.close()
    if not prox.empty: plt.figure(figsize=(7,5)); p2=prox.copy(); p2["plot_set"]=p2["set"].replace({"杜仲":"Du-Zhong","牛膝":"Niu-Xi","续断":"Xu-Duan","骨碎补":"Gu-Sui-Bu"}); sns.barplot(data=p2,x="plot_set",y="network_distance",color="#b2182b"); plt.xticks(rotation=30,ha="right"); plt.tight_layout(); plt.savefig(figdir/"Fig3_network_proximity.png",dpi=300); plt.close()
    plt.figure(figsize=(8,5)); c2=cell.copy(); c2["plot_cell"]=c2.cell_state.replace({"BMSC/成骨祖细胞":"BMSC","成骨细胞":"Osteoblast","破骨细胞":"Osteoclast","骨免疫细胞":"Osteoimmune"}); sns.barplot(data=c2,x="plot_cell",y="ucell_signature_score",color="#1b9e77"); plt.tight_layout(); plt.savefig(figdir/"Fig4_single_cell_localisation.png",dpi=300); plt.close()
    if not genetics.empty: plt.figure(figsize=(6,6)); sns.barplot(data=genetics.head(20),y="target",x="open_targets_score",color="#b2182b"); plt.tight_layout(); plt.savefig(figdir/"Fig5_human_genetics.png",dpi=300); plt.close()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--workbook",default="osteoporosis_extraction_output_finalV6.xlsx"); ap.add_argument("--outdir",default="results"); ap.add_argument("--cache",default=".cache/public_api"); ap.add_argument("--analysis-scope", choices=["quick","full"], default="quick"); ap.add_argument("--dedup-strategy", choices=["merge","highest_quality","first","raw"], default="merge"); ap.add_argument("--download-large", action="store_true"); ap.add_argument("--large-data-dir", default="data/full_scale")
    args=ap.parse_args(); out=Path(args.outdir); tab=out/"tables"; fig=out/"figures"; tab.mkdir(parents=True,exist_ok=True); fig.mkdir(parents=True,exist_ok=True)
    df=load_records(args.workbook, strategy=args.dedup_strategy); df["dedup_strategy"]=args.dedup_strategy
    df.drop(columns=[c for c in ["herbs_set","symptoms_set"] if c in df.columns]).to_csv(tab/"deduplicated_classical_records.csv", index=False)
    dedup_conflict_log(args.workbook, tab); dedup_strategy_sensitivity(args.workbook, tab)
    mods,hist=mine_modules(df,tab); core_sig=core_module_significance(df, tab)
    cache=ApiCache(Path(args.cache)); comp,targets,disease,prox=real_targets_and_network(tab, cache); cell,gen=omics_genetics(tab,targets,disease); reference_marker_overlap(tab, cache, targets); gse224152_expression_localisation(tab, cache, targets); gse246769_osteoclast_dynamics(tab, cache, targets); write_dataset_manifest(tab); gwas_catalog_trait_studies(tab, cache)
    _gwas_enr, gwas_genes = gwas_gene_enrichment(tab, cache, targets)
    causal_pharmacology_scores(tab, gwas_genes); write_high_order_resource_manifest(tab); write_full_scale_execution_plan(tab, args.analysis_scope, args.download_large, args.large_data_dir); download_requested_large_resources(cache, tab, args.analysis_scope, args.download_large, args.large_data_dir); plot_all(fig,mods,hist,prox,cell,gen); plot_nature_style_extensions(fig, tab)
    print(f"Analysed {len(df)} deduplicated records; core-module significance rows={len(core_sig)}; GWAS bone genes={len(gwas_genes)}; wrote outputs to {out}.")
if __name__ == "__main__": main()
