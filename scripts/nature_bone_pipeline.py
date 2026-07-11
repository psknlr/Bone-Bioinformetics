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
from scipy.stats import fisher_exact
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

def load_records(path):
    df = pd.read_excel(path, sheet_name="Sheet1")
    df = df[df["status"].eq("ok") & df["consider_include"].fillna(False)].copy()
    df["record_key"] = df[["Books","year","Diagnosis","Chapter","Title"]].fillna("").astype(str).agg("|".join, axis=1).map(lambda s: hashlib.md5(s.encode()).hexdigest())
    df = df.drop_duplicates("record_key")
    df["herbs"] = [[h for h in sorted({norm_herb(d.get("name", "")) for d in parse_json_list(x) if d.get("name")}) if h] for x in df["therapeutic_drugs_json"]]
    df = df[df["herbs"].map(len).gt(0)].copy()
    df["dynasty"] = pd.cut(df["year"], bins=[x[0] for x in DYNASTIES]+[10_000], labels=[x[2] for x in DYNASTIES], right=False)
    return df

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
    for group_col in ["dynasty","tcm_syndrome","Diagnosis"]:
        for g, sub in df.groupby(group_col, dropna=True):
            hist.append({"stratum":group_col,"level":str(g),"records":len(sub),"core_complete_n":sum(core.issubset(set(x)) for x in sub.herbs),"core_any_n":sum(bool(core & set(x)) for x in sub.herbs),"core_any_rate":sum(bool(core & set(x)) for x in sub.herbs)/len(sub)})
    pd.DataFrame(hist).to_csv(outdir/"historical_stability.csv", index=False); return mods, pd.DataFrame(hist)

def pubchem_lookup(cache, name):
    url=f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/property/Title,MolecularFormula,CanonicalSMILES/JSON"
    try: props=cache.get_json(url)["data"]["PropertyTable"]["Properties"][0]
    except Exception: return None
    return {"query":name,"pubchem_cid":props.get("CID"),"pubchem_title":props.get("Title"),"formula":props.get("MolecularFormula"),"canonical_smiles":props.get("CanonicalSMILES")}

def chembl_molecule(cache, name):
    data=cache.get_json("https://www.ebi.ac.uk/chembl/api/data/molecule/search.json", {"q":name,"limit":1})["data"]
    mols=data.get("molecules") or []
    return mols[0].get("molecule_chembl_id") if mols else None

def chembl_targets(cache, chembl_id, min_pchembl=5.0, max_pages=4):
    out=[]
    for page in range(max_pages):
        data=cache.get_json("https://www.ebi.ac.uk/chembl/api/data/activity.json", {"molecule_chembl_id":chembl_id,"pchembl_value__isnull":False,"limit":100,"offset":page*100})["data"]
        for a in data.get("activities", []):
            try: pchem=float(a.get("pchembl_value") or 0)
            except ValueError: pchem=0
            if pchem >= min_pchembl and a.get("target_organism") == "Homo sapiens":
                out.append({"molecule_chembl_id":chembl_id,"target_chembl_id":a.get("target_chembl_id"),"target_pref_name":a.get("target_pref_name"),"pchembl_value":pchem,"standard_type":a.get("standard_type")})
        if not data.get("page_meta",{}).get("next"): break
    return out

def chembl_target_gene_symbol(cache, target_chembl_id):
    if not target_chembl_id: return None
    try:
        data=cache.get_json(f"https://www.ebi.ac.uk/chembl/api/data/target/{target_chembl_id}.json")["data"]
    except Exception:
        return None
    for comp in data.get("target_components") or []:
        for syn in comp.get("target_component_synonyms") or []:
            if syn.get("syn_type") == "GENE_SYMBOL":
                return syn.get("component_synonym")
    return None

def opentargets_osteoporosis(cache, size=200):
    q='''query disease($id:String!,$size:Int!){ disease(efoId:$id){ id name associatedTargets(page:{index:0,size:$size}){ rows{ score target{ approvedSymbol id } } } } }'''
    data=cache.get_json("https://api.platform.opentargets.org/api/v4/graphql", method="POST", payload={"query":q,"variables":{"id":"MONDO_0005298","size":size}})["data"]
    rows=data.get("data",{}).get("disease",{}).get("associatedTargets",{}).get("rows",[]) if data else []
    return pd.DataFrame([{"target":r["target"].get("approvedSymbol"),"open_targets_score":r.get("score"),"ensembl_id":r["target"].get("id")} for r in rows if r.get("target")])

def string_network(cache, genes):
    genes=sorted({g for g in genes if isinstance(g,str) and g});
    if len(genes)<2: return nx.Graph()
    text=requests.get("https://string-db.org/api/tsv/network", params={"identifiers":"%0d".join(genes),"species":9606,"required_score":400}, timeout=60).text
    G=nx.Graph(); G.add_nodes_from(genes)
    for line in text.splitlines()[1:]:
        f=line.split('\t')
        if len(f)>5: G.add_edge(f[2], f[3], score=float(f[5]))
    return G

def real_targets_and_network(outdir, cache):
    compounds=[]; target_rows=[]
    for herb, queries in HERB_COMPOUND_QUERIES.items():
        for q in queries:
            pc=pubchem_lookup(cache, q) or {"query":q,"pubchem_cid":None,"pubchem_title":None,"formula":None,"canonical_smiles":None}
            chembl=chembl_molecule(cache, q); pc.update({"herb":herb,"molecule_chembl_id":chembl}); compounds.append(pc)
            if chembl:
                for t in chembl_targets(cache, chembl):
                    t.update({"herb":herb,"compound_query":q,"target_gene_symbol":chembl_target_gene_symbol(cache, t.get("target_chembl_id"))}); target_rows.append(t)
    comp=pd.DataFrame(compounds); comp.to_csv(outdir/"pubchem_chembl_compounds.csv", index=False)
    targets=pd.DataFrame(target_rows); targets.to_csv(outdir/"component_target_evidence_tiers.csv", index=False)
    disease=opentargets_osteoporosis(cache); disease.to_csv(outdir/"opentargets_osteoporosis_targets.csv", index=False)
    drug_gene=set(targets.get("target_gene_symbol", pd.Series(dtype=str)).dropna())
    disease_gene=set(disease.target.dropna().head(100)) if not disease.empty else set()
    G=string_network(cache, drug_gene | disease_gene)
    def prox(genes):
        vals=[]
        for g in genes:
            if g in G:
                ds=[nx.shortest_path_length(G,g,d) for d in disease_gene if d in G and nx.has_path(G,g,d) and g!=d]
                if ds: vals.append(min(ds))
        return np.mean(vals) if vals else np.nan
    rows=[{"set":"core_combo","targets":len(drug_gene),"network_distance":prox(drug_gene),"source":"ChEMBL+OpenTargets+STRING"}]
    for herb in CORE:
        hg=set(targets.loc[targets.herb.eq(herb),"target_gene_symbol"].dropna()) if not targets.empty else set()
        rows.append({"set":herb,"targets":len(hg),"network_distance":prox(hg),"source":"ChEMBL+OpenTargets+STRING"})
    pd.DataFrame(rows).to_csv(outdir/"network_proximity.csv", index=False)
    return comp, targets, disease, pd.DataFrame(rows)

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
    # PanglaoDB: real downloaded marker table.
    try:
        raw=cache.get_bytes(PANGLOADB_URL, headers={"User-Agent":"Mozilla/5.0"})
        pang=pd.read_csv(io.BytesIO(gzip.decompress(raw)), sep="\t")
        keep=pang["cell type"].str.lower().isin({"osteoblasts","osteoclasts","osteoclast precursor cells","osteocytes","stromal cells"})
        pang=pang[keep & pang["species"].str.contains("Hs", na=False)].copy()
        for cell, sub in pang.groupby("cell type"):
            markers=set(sub["official gene symbol"].dropna().astype(str))
            marker_rows.append({"source":"PanglaoDB_27_Mar_2020","cell_state":cell,"marker_n":len(markers),"target_overlap":len(genes & markers),"overlap_genes":";".join(sorted(genes & markers))})
    except Exception as e:
        marker_rows.append({"source":"PanglaoDB_27_Mar_2020","cell_state":"download_failed","marker_n":0,"target_overlap":0,"overlap_genes":str(e)[:200]})
    for cell, markers in CELL_MARKERS.items():
        marker_rows.append({"source":"curated_bone_marker_panel","cell_state":cell,"marker_n":len(markers),"target_overlap":len(genes & markers),"overlap_genes":";".join(sorted(genes & markers))})
    out=pd.DataFrame(marker_rows); out.to_csv(outdir/"reference_marker_overlap_localisation.csv", index=False)
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

def causal_pharmacology_scores(outdir):
    targets=pd.read_csv(outdir/"component_target_evidence_tiers.csv")
    ot=pd.read_csv(outdir/"opentargets_osteoporosis_targets.csv") if (outdir/"opentargets_osteoporosis_targets.csv").exists() else pd.DataFrame(columns=["target","open_targets_score"])
    marker=pd.read_csv(outdir/"reference_marker_overlap_localisation.csv") if (outdir/"reference_marker_overlap_localisation.csv").exists() else pd.DataFrame(columns=["overlap_genes"])
    gse224=pd.read_csv(outdir/"gse224152_target_expression.csv") if (outdir/"gse224152_target_expression.csv").exists() else pd.DataFrame(columns=["gene","pct_cells_detected","mean_expression"])
    gse246=pd.read_csv(outdir/"gse246769_osteoclast_dynamics.csv") if (outdir/"gse246769_osteoclast_dynamics.csv").exists() else pd.DataFrame(columns=["gene","delta_d9_vs_d0"])
    marker_genes=set()
    for x in marker.get("overlap_genes", pd.Series(dtype=str)).dropna():
        marker_genes |= {g for g in str(x).split(";") if g and g != "nan"}
    rows=[]
    for gene, sub in targets.dropna(subset=["target_gene_symbol"]).groupby("target_gene_symbol"):
        ot_score=float(ot.loc[ot.target.eq(gene),"open_targets_score"].max()) if gene in set(ot.target) else 0.0
        pct=float(gse224.loc[gse224.gene.eq(gene),"pct_cells_detected"].max()) if gene in set(gse224.gene) else 0.0
        delta=float(gse246.loc[gse246.gene.eq(gene),"delta_d9_vs_d0"].mean()) if gene in set(gse246.gene) else 0.0
        max_pchem=float(sub.pchembl_value.max()) if "pchembl_value" in sub else 0.0
        evidence_count=sum([ot_score>0, gene in marker_genes, pct>0.01, abs(delta)>1, max_pchem>=6])
        if ot_score>0 and evidence_count>=4: grade="A_candidate_requires_coloc_MR_confirmation"
        elif ot_score>0 or evidence_count>=3: grade="B_multi_omics_partial_support"
        else: grade="C_pharmacology_only_or_weak_omics"
        direction="inhibit_candidate" if delta>1 else ("activate_candidate" if delta<-1 or pct>0.05 else "direction_uncertain")
        rows.append({"target":gene,"herbs":";".join(sorted(sub.herb.dropna().unique())),"max_pchembl":max_pchem,"open_targets_score":ot_score,"reference_marker_overlap":gene in marker_genes,"gse224152_pct_cells_detected":pct,"gse246769_delta_d9_vs_d0":delta,"evidence_count":evidence_count,"causal_evidence_grade":grade,"estimated_pharmacology_direction":direction,"fine_mapping_coloc_mr_status":"not_run_requires_full_sumstats_QTL"})
    out=pd.DataFrame(rows).sort_values(["causal_evidence_grade","open_targets_score","evidence_count","max_pchembl"], ascending=[True,False,False,False])
    out.to_csv(outdir/"genetics_anchored_causal_pharmacology_scores.csv", index=False); return out

def write_high_order_resource_manifest(outdir):
    out=pd.DataFrame(HIGH_ORDER_RESOURCES); out.to_csv(outdir/"high_order_public_resource_manifest.csv", index=False); return out

def plot_all(figdir, mods, hist, prox, cell, genetics):
    sns.set_theme(style="whitegrid", context="talk")
    plt.figure(figsize=(10,7)); top=mods.nsmallest(20,"fdr").copy(); top["-log10(FDR)"]=-np.log10(top.fdr.clip(1e-12)); sns.scatterplot(data=top,x="lift_vs_independence",y="-log10(FDR)",size="support_n",hue="contains_core",sizes=(80,450)); plt.tight_layout(); plt.savefig(figdir/"Fig1_stable_modules.png",dpi=300); plt.close()
    plt.figure(figsize=(10,6)); h=hist[hist.stratum.eq("dynasty")].copy(); h["plot_level"]=h.level.map({"隋唐以前":"Pre-Tang","唐":"Tang","宋金元":"Song-Jin-Yuan","明":"Ming","清":"Qing","近现代":"Modern"}); sns.barplot(data=h,x="plot_level",y="core_any_rate",color="#2166ac"); plt.xticks(rotation=35,ha="right"); plt.tight_layout(); plt.savefig(figdir/"Fig2_historical_stability.png",dpi=300); plt.close()
    if not prox.empty: plt.figure(figsize=(7,5)); p2=prox.copy(); p2["plot_set"]=p2["set"].replace({"杜仲":"Du-Zhong","牛膝":"Niu-Xi","续断":"Xu-Duan","骨碎补":"Gu-Sui-Bu"}); sns.barplot(data=p2,x="plot_set",y="network_distance",color="#b2182b"); plt.xticks(rotation=30,ha="right"); plt.tight_layout(); plt.savefig(figdir/"Fig3_network_proximity.png",dpi=300); plt.close()
    plt.figure(figsize=(8,5)); c2=cell.copy(); c2["plot_cell"]=c2.cell_state.replace({"BMSC/成骨祖细胞":"BMSC","成骨细胞":"Osteoblast","破骨细胞":"Osteoclast","骨免疫细胞":"Osteoimmune"}); sns.barplot(data=c2,x="plot_cell",y="ucell_signature_score",color="#1b9e77"); plt.tight_layout(); plt.savefig(figdir/"Fig4_single_cell_localisation.png",dpi=300); plt.close()
    if not genetics.empty: plt.figure(figsize=(6,6)); sns.barplot(data=genetics.head(20),y="target",x="open_targets_score",color="#b2182b"); plt.tight_layout(); plt.savefig(figdir/"Fig5_human_genetics.png",dpi=300); plt.close()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--workbook",default="osteoporosis_extraction_output_finalV6.xlsx"); ap.add_argument("--outdir",default="results"); ap.add_argument("--cache",default=".cache/public_api")
    args=ap.parse_args(); out=Path(args.outdir); tab=out/"tables"; fig=out/"figures"; tab.mkdir(parents=True,exist_ok=True); fig.mkdir(parents=True,exist_ok=True)
    df=load_records(args.workbook); df.to_csv(tab/"deduplicated_classical_records.csv", index=False)
    mods,hist=mine_modules(df,tab); cache=ApiCache(Path(args.cache)); comp,targets,disease,prox=real_targets_and_network(tab, cache); cell,gen=omics_genetics(tab,targets,disease); reference_marker_overlap(tab, cache, targets); gse224152_expression_localisation(tab, cache, targets); gse246769_osteoclast_dynamics(tab, cache, targets); write_dataset_manifest(tab); gwas_catalog_trait_studies(tab, cache); causal_pharmacology_scores(tab); write_high_order_resource_manifest(tab); plot_all(fig,mods,hist,prox,cell,gen)
    print(f"Analysed {len(df)} deduplicated records; queried public APIs; wrote outputs to {out}.")
if __name__ == "__main__": main()
