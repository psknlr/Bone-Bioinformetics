#!/usr/bin/env python3
"""High-order public-data modules that extend the core osteoporosis herb pipeline.

Design contract (identical to scripts/nature_bone_pipeline.py):
  * every step that a public REST/GraphQL endpoint supports without multi-GB
    downloads is computed for real and cached;
  * steps that genuinely require full GWAS summary statistics, QTL tabix files,
    large single-cell/spatial matrices, or docking/MD binaries are written as
    explicit, structured `pending` records naming the exact required input —
    they are never simulated.

Modules
  1. Genetics-anchored multi-omics causal pharmacology (Open Targets genetic
     datatype evidence + GWAS-Catalog gene overlap + eQTL presence; fine-mapping/
     coloc/SMR/MR marked pending on full sumstats/QTL).
  2. Single-cell / spatial bone-niche localisation (real module-expression and
     marker-niche assignment on downloaded matrices; scVI/cell2location/CellChat/
     trajectory pending on large matrices).
  3. Perturbation-omics disease-signature reversal (real Enrichr LINCS L1000
     query from a GSE246769-derived osteoclast signature; Level-5 weighted CMap
     and scPerturb E-distance pending on bulk signature files).
  4. Heterogeneous knowledge graph + structural pharmacology (real KG from
     pipeline tables + Reactome pathways + Open Targets, meta-path link
     prediction, UniProt/AlphaFold/PDB structure availability; docking/MD pending).

Run after nature_bone_pipeline.py has written results/tables.
"""
from __future__ import annotations

import argparse, io, time
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
import requests

from nature_bone_pipeline import ApiCache, CORE

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
OSTEOPOROSIS_EFO = "MONDO_0005298"

# Four spatial pharmacology niches named in the study design, with marker panels.
NICHE_MARKERS = {
    "osteoprogenitor_niche": {"LEPR", "CXCL12", "PDGFRA", "RUNX2", "SP7", "BMP2", "WNT16", "LRP5", "COL1A1"},
    "osteoclast_immune_niche": {"TNFSF11", "TNFRSF11A", "CTSK", "ACP5", "NFATC1", "SPP1", "MMP9", "CSF1R"},
    "vascular_osteogenic_niche": {"PECAM1", "CDH5", "KDR", "VEGFA", "ANGPT1", "NOTCH1", "ESM1"},
    "marrow_adipogenic_niche": {"PPARG", "CEBPA", "ADIPOQ", "LPL", "FABP4", "CFD"},
}


# --------------------------------------------------------------------------- #
# Module 1 — genetics-anchored multi-omics causal pharmacology
# --------------------------------------------------------------------------- #
def opentargets_genetic_datatypes(cache, size=400):
    """Open Targets datatype scores per osteoporosis-associated target; the
    `genetic_association` component is the real human-genetics anchor."""
    q = ("query d($id:String!,$size:Int!){disease(efoId:$id){associatedTargets(page:{index:0,size:$size})"
         "{rows{score target{approvedSymbol id} datatypeScores{id score}}}}}")
    data = cache.get_json(OT_GRAPHQL, method="POST",
                          payload={"query": q, "variables": {"id": OSTEOPOROSIS_EFO, "size": size}})["data"]
    rows = (data.get("data", {}).get("disease", {}) or {}).get("associatedTargets", {}).get("rows", []) if data else []
    out = []
    for r in rows:
        tgt = r.get("target") or {}
        dts = {d["id"]: d["score"] for d in r.get("datatypeScores", [])}
        out.append({"target": tgt.get("approvedSymbol"), "ensembl_id": tgt.get("id"),
                    "overall_association_score": r.get("score"),
                    "genetic_association": dts.get("genetic_association", 0.0),
                    "genetic_literature": dts.get("genetic_literature", 0.0),
                    "rna_expression": dts.get("rna_expression", 0.0),
                    "animal_model": dts.get("animal_model", 0.0),
                    "known_drug": dts.get("known_drug", 0.0)})
    return pd.DataFrame(out)


def symbol_to_ensembl(cache, symbols):
    """Map gene symbols to Ensembl IDs through MyGene.info (real, per-symbol GET, cached)."""
    mapping = {}
    for s in sorted({s for s in symbols if isinstance(s, str) and s}):
        try:
            data = cache.get_json("https://mygene.info/v3/query",
                                  {"q": f"symbol:{s}", "species": "human", "fields": "ensembl.gene", "size": 1})["data"]
        except Exception:
            continue
        hits = data.get("hits", []) if isinstance(data, dict) else []
        if hits:
            ens = hits[0].get("ensembl")
            gid = ens.get("gene") if isinstance(ens, dict) else (ens[0].get("gene") if isinstance(ens, list) and ens else None)
            if gid:
                mapping[s] = gid
    return mapping


def _eqtl_bone_relevant_datasets(cache, want=6):
    """A small curated set of eQTL Catalogue gene-expression datasets in
    bone/marrow-relevant tissues (blood, macrophage, monocyte, muscle, fibroblast)."""
    try:
        dsets = cache.get_json("https://www.ebi.ac.uk/eqtl/api/v2/datasets", {"size": 1000})["data"]
    except Exception:
        return []
    keep_terms = ("blood", "macrophage", "monocyte", "muscle", "fibroblast", "mesenchym", "osteo")
    out = []
    for d in dsets if isinstance(dsets, list) else []:
        if d.get("quant_method") != "ge":
            continue
        if any(t in str(d.get("tissue_label", "")).lower() for t in keep_terms):
            out.append(d["dataset_id"])
        if len(out) >= want:
            break
    return out

def eqtl_catalogue_cis_support(cache, ensembl_ids, max_genes=60):
    """Real cis-eQTL presence check: for each gene, query a curated set of
    bone/immune-relevant eQTL Catalogue datasets (per-dataset REST path) and mark
    True if any significant (p<5e-8) cis-eQTL is found, False if none, None if the
    catalogue is unreachable. Full multi-tissue coloc/SMR/MR stay pending."""
    datasets = _eqtl_bone_relevant_datasets(cache)
    rows = []
    for gid in list(ensembl_ids)[:max_genes]:
        if not datasets:
            rows.append({"ensembl_id": gid, "has_significant_cis_eqtl": None, "n_datasets_queried": 0}); continue
        has = False
        for ds in datasets:
            try:
                data = cache.get_json(f"https://www.ebi.ac.uk/eqtl/api/v2/datasets/{ds}/associations",
                                      {"gene_id": gid, "size": 1, "p_upper": 5e-8})["data"]
                if (data if isinstance(data, list) else data.get("_embedded", {}).get("associations", [])):
                    has = True; break
            except Exception:
                continue  # 400 "No results" or transient -> no hit in this dataset
        rows.append({"ensembl_id": gid, "has_significant_cis_eqtl": has, "n_datasets_queried": len(datasets)})
    return pd.DataFrame(rows)


def module1_causal_pharmacology(outdir, cache):
    tiers = pd.read_csv(outdir / "component_target_evidence_tiers.csv")
    if "target_gene_symbol" not in tiers.columns:
        return
    genes = sorted({str(g) for g in tiers["target_gene_symbol"].dropna()})
    gwas_bone = set()
    p = outdir / "gwas_catalog_bone_genes.csv"
    if p.exists():
        gwas_bone = {str(g).upper() for g in pd.read_csv(p)["gene"].dropna()}
    ot = opentargets_genetic_datatypes(cache)
    ot_gen = {r.target: r for _, r in ot.iterrows()} if not ot.empty else {}
    ens_map = symbol_to_ensembl(cache, genes)
    eqtl = eqtl_catalogue_cis_support(cache, set(ens_map.values()))
    eqtl_by_ens = {r.ensembl_id: r.has_significant_cis_eqtl for _, r in eqtl.iterrows()}
    # direction proxy from real osteoclast dynamics
    dyn = pd.read_csv(outdir / "gse246769_osteoclast_dynamics.csv") if (outdir / "gse246769_osteoclast_dynamics.csv").exists() else pd.DataFrame(columns=["gene", "delta_d9_vs_d0"])
    rows = []
    for gene, sub in tiers.dropna(subset=["target_gene_symbol"]).groupby("target_gene_symbol"):
        rec = ot_gen.get(gene)
        gen_score = float(rec["genetic_association"]) if rec is not None else 0.0
        in_gwas = gene.upper() in gwas_bone
        ens = ens_map.get(gene)
        has_eqtl = eqtl_by_ens.get(ens)
        max_pchem = float(sub["pchembl_value"].max()) if "pchembl_value" in sub else 0.0
        delta = float(dyn.loc[dyn.gene.eq(gene), "delta_d9_vs_d0"].mean()) if gene in set(dyn.gene) else np.nan
        # genetics-anchored grade: A requires real human-genetics support
        has_human_genetics = (gen_score >= 0.10) or in_gwas
        strong_genetics = (gen_score >= 0.10) and (in_gwas or has_eqtl is True)
        if strong_genetics and max_pchem >= 6:
            grade = "A_genetic_finemap_coloc_MR_priority"
        elif has_human_genetics:
            grade = "B_partial_genetics_or_omics_support"
        else:
            grade = "C_pharmacology_or_prediction_only"
        # osteoclast-differentiation direction → anti-resorptive vs anabolic hypothesis
        if not np.isnan(delta) and delta > 1:
            direction = "up_in_osteoclastogenesis:inhibit_to_reduce_resorption"
        elif not np.isnan(delta) and delta < -1:
            direction = "down_in_osteoclastogenesis:activate_or_context_dependent"
        else:
            direction = "direction_uncertain_requires_MR_beta"
        rows.append({"target": gene, "herbs": ";".join(sorted(sub.herb.dropna().unique())),
                     "ot_genetic_association": gen_score, "gwas_catalog_bone_gene": in_gwas,
                     "has_significant_cis_eqtl": has_eqtl, "max_pchembl": max_pchem,
                     "osteoclast_delta_d9_vs_d0": delta, "human_genetics_support": has_human_genetics,
                     "causal_grade": grade, "pharmacology_direction_hypothesis": direction})
    res = pd.DataFrame(rows).sort_values(["causal_grade", "ot_genetic_association", "max_pchembl"],
                                         ascending=[True, False, False])
    res.to_csv(outdir / "m1_genetics_anchored_causal_targets.csv", index=False)
    pending = pd.DataFrame([
        {"step": "SuSiE/FINEMAP fine-mapping", "requires": "full GWAS summary statistics (GWAS Catalog / OT credibleSets)", "status": "pending_inputs"},
        {"step": "coloc / SuSiE-coloc", "requires": "GWAS region sumstats + eQTL Catalogue tabix", "status": "pending_inputs"},
        {"step": "SMR/HEIDI", "requires": "GWAS sumstats + cis-eQTL/pQTL besd", "status": "pending_inputs"},
        {"step": "cis-eQTL/pQTL Mendelian randomisation", "requires": "harmonised instrument sumstats", "status": "pending_inputs"},
        {"step": "multi-tissue QTL consistency", "requires": "eQTL Catalogue multi-dataset tabix", "status": "pending_inputs"},
    ])
    pending.to_csv(outdir / "m1_causal_pending_steps.csv", index=False)
    ot.to_csv(outdir / "m1_opentargets_genetic_datatypes.csv", index=False)
    return res


# --------------------------------------------------------------------------- #
# Module 2 — single-cell / spatial bone-niche localisation
# --------------------------------------------------------------------------- #
def module2_niche_localisation(outdir):
    tiers = pd.read_csv(outdir / "component_target_evidence_tiers.csv")
    genes = {str(g) for g in tiers.get("target_gene_symbol", pd.Series(dtype=str)).dropna()}
    # marker-niche assignment (real overlap with named niches)
    rows = []
    for niche, markers in NICHE_MARKERS.items():
        ov = genes & markers
        rows.append({"niche": niche, "marker_n": len(markers), "target_overlap": len(ov),
                     "overlap_genes": ";".join(sorted(ov)),
                     "assignment": "candidate_niche" if ov else "no_direct_marker_overlap"})
    pd.DataFrame(rows).to_csv(outdir / "m2_niche_localisation.csv", index=False)
    # real module-expression scores on the downloaded osteoclast time series
    dyn_p = outdir / "gse246769_osteoclast_dynamics.csv"
    if dyn_p.exists():
        dyn = pd.read_csv(dyn_p)
        cols = [c for c in ["log2cpm_d0", "log2cpm_d2", "log2cpm_d5", "log2cpm_d9"] if c in dyn.columns]
        mod = pd.DataFrame([{"dataset": "GSE246769", "day": c.replace("log2cpm_", ""),
                             "module_mean_log2cpm": float(dyn[c].mean()),
                             "n_targets_scored": int(dyn[c].notna().sum())} for c in cols])
        mod.to_csv(outdir / "m2_module_expression_scores.csv", index=False)
    pending = pd.DataFrame([
        {"step": "scVI/scANVI or Harmony cross-dataset integration", "requires": "GSE169396 + GSE253355 + GSE284089 matrices", "status": "pending_inputs"},
        {"step": "donor-level pseudo-bulk differential", "requires": "annotated single-cell matrices with donor labels", "status": "pending_inputs"},
        {"step": "cell2location / RCTD / Tangram spatial mapping", "requires": "GSE284089 spatial slides + single-cell reference", "status": "pending_inputs"},
        {"step": "UCell / AUCell single-cell module score", "requires": "single-cell expression matrix", "status": "pending_inputs"},
        {"step": "Monocle3 / Slingshot / CellRank trajectory", "requires": "BMSC→osteoblast / monocyte→osteoclast matrices", "status": "pending_inputs"},
        {"step": "LIANA+/CellChat/NicheNet communication", "requires": "annotated single-cell matrices", "status": "pending_inputs"},
        {"step": "pySCENIC / DoRothEA regulon inference", "requires": "single-cell counts", "status": "pending_inputs"},
    ])
    pending.to_csv(outdir / "m2_pending_steps.csv", index=False)


# --------------------------------------------------------------------------- #
# Module 3 — perturbation-omics disease-signature reversal (Enrichr LINCS)
# --------------------------------------------------------------------------- #
def enrichr_query(cache, genes, library):
    genes = [g for g in genes if isinstance(g, str) and g]
    if len(genes) < 3:
        return pd.DataFrame()
    r = requests.post("https://maayanlab.cloud/Enrichr/addList",
                      files={"list": (None, "\n".join(genes)), "description": (None, "sig")}, timeout=45)
    r.raise_for_status(); uid = r.json()["userListId"]; time.sleep(0.3)
    r2 = requests.get("https://maayanlab.cloud/Enrichr/enrich",
                      params={"userListId": uid, "backgroundType": library}, timeout=60)
    r2.raise_for_status(); terms = r2.json().get(library, [])
    return pd.DataFrame([{"library": library, "term": t[1], "p_value": t[2], "combined_score": t[4],
                          "overlap_genes": ";".join(t[5])} for t in terms[:100]])


def module3_perturbation_reversal(outdir, cache):
    dyn_p = outdir / "gse246769_osteoclast_dynamics.csv"
    pending = pd.DataFrame([
        {"step": "weighted L1000 Level-5 connectivity score", "requires": "GSE92742/GSE70138 Level-5 signatures", "status": "pending_inputs"},
        {"step": "scPerturb E-distance effect comparison", "requires": "scPerturb h5ad datasets", "status": "pending_inputs"},
        {"step": "combination Pareto-front synergy modelling", "requires": "per-compound reversal profiles", "status": "pending_inputs"},
    ])
    pending.to_csv(outdir / "m3_pending_steps.csv", index=False)
    if not dyn_p.exists():
        return
    dyn = pd.read_csv(dyn_p)
    up = sorted(dyn.loc[dyn.delta_d9_vs_d0 > 1, "gene"].dropna().unique())
    down = sorted(dyn.loc[dyn.delta_d9_vs_d0 < -1, "gene"].dropna().unique())
    frames = []
    for lib in ["LINCS_L1000_Chem_Pert_Consensus_Sigs", "LINCS_L1000_CRISPR_KO_Consensus_Sigs"]:
        for sig_name, sig in [("osteoclast_up_at_d9", up), ("osteoclast_down_at_d9", down)]:
            try:
                res = enrichr_query(cache, sig, lib)
            except Exception as e:
                res = pd.DataFrame([{"library": lib, "term": f"query_failed:{str(e)[:80]}", "p_value": np.nan, "combined_score": np.nan, "overlap_genes": ""}])
            if not res.empty:
                res.insert(0, "disease_signature", sig_name); frames.append(res)
    if frames:
        out = pd.concat(frames, ignore_index=True)
        out["interpretation"] = "perturbagens associated with the osteoclast signature; candidate reversers where a perturbagen down-regulates the up-signature (confirm with Level-5 directionality)"
        out.to_csv(outdir / "m3_perturbation_reversal.csv", index=False)


# --------------------------------------------------------------------------- #
# Module 4 — heterogeneous knowledge graph + structural pharmacology
# --------------------------------------------------------------------------- #
def reactome_pathways(cache, uniprot):
    try:
        url = f"https://reactome.org/ContentService/data/mapping/UniProt/{uniprot}/pathways"
        data = cache.get_json(url, {"species": "9606"})["data"]
        return [(p.get("stId"), p.get("displayName")) for p in data] if isinstance(data, list) else []
    except Exception:
        return []


def uniprot_accession(cache, symbol):
    try:
        data = cache.get_json("https://rest.uniprot.org/uniprotkb/search",
                              {"query": f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true",
                               "fields": "accession", "size": 1, "format": "json"})["data"]
        res = data.get("results", [])
        return res[0].get("primaryAccession") if res else None
    except Exception:
        return None


def alphafold_available(cache, uniprot):
    try:
        data = cache.get_json(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}")["data"]
        if data:
            return True, data[0].get("globalMetricValue")
    except Exception:
        pass
    return None, None


def pdb_available(cache, uniprot):
    try:
        data = cache.get_json(f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot}")["data"]
        recs = data.get(uniprot, []) if isinstance(data, dict) else []
        return bool(recs), (recs[0].get("pdb_id") if recs else None)
    except Exception:
        return None, None


def module4_knowledge_graph(outdir, cache, top_structure=20):
    tiers = pd.read_csv(outdir / "component_target_evidence_tiers.csv")
    comp = pd.read_csv(outdir / "pubchem_chembl_compounds.csv") if (outdir / "pubchem_chembl_compounds.csv").exists() else pd.DataFrame()
    ot = pd.read_csv(outdir / "opentargets_osteoporosis_targets.csv") if (outdir / "opentargets_osteoporosis_targets.csv").exists() else pd.DataFrame(columns=["target"])
    disease_genes = set(ot.target.dropna()) if not ot.empty else set()
    G = nx.Graph(); edges = []
    def add(a, ta, b, tb, rel):
        G.add_node(a, ntype=ta); G.add_node(b, ntype=tb); G.add_edge(a, b, rel=rel)
        edges.append({"source": a, "source_type": ta, "target": b, "target_type": tb, "relation": rel})
    for _, r in tiers.dropna(subset=["target_gene_symbol"]).iterrows():
        herb, cq, gene = r.get("herb"), r.get("compound_query"), r["target_gene_symbol"]
        if herb and cq: add(str(herb), "herb", str(cq), "compound", "contains")
        if cq and gene: add(str(cq), "compound", str(gene), "target", "binds")
    for gene in disease_genes:
        add(str(gene), "target", "osteoporosis", "disease", "associated_with")
    # Reactome pathway layer for candidate targets
    genes = sorted({str(g) for g in tiers.get("target_gene_symbol", pd.Series(dtype=str)).dropna()})
    uni = {}
    for g in genes[:top_structure * 3]:
        acc = uniprot_accession(cache, g)
        if acc:
            uni[g] = acc
            for stid, name in reactome_pathways(cache, acc)[:8]:
                add(str(g), "target", str(name), "pathway", "participates_in")
    pd.DataFrame(edges).to_csv(outdir / "m4_knowledge_graph_edges.csv", index=False)
    # meta-path link prediction: rank herb→disease-gene candidates by shared-pathway / Adamic-Adar
    target_nodes = [n for n, d in G.nodes(data=True) if d.get("ntype") == "target"]
    preds = []
    for herb in CORE:
        if herb not in G:
            continue
        herb_targets = {nb for nb in nx.single_source_shortest_path_length(G, herb, cutoff=2) if G.nodes.get(nb, {}).get("ntype") == "target"}
        for dg in disease_genes:
            if dg not in G or dg in herb_targets:
                continue
            # shared pathway neighbours between herb's targets and the disease gene
            dg_paths = {nb for nb in G.neighbors(dg) if G.nodes[nb].get("ntype") == "pathway"}
            shared = 0
            for ht in herb_targets:
                shared += len(dg_paths & {nb for nb in G.neighbors(ht) if G.nodes[nb].get("ntype") == "pathway"})
            if shared:
                preds.append({"herb": herb, "candidate_target": dg, "shared_pathway_links": shared,
                              "via_experimental_targets": ";".join(sorted(herb_targets & set(target_nodes))[:8])})
    pd.DataFrame(sorted(preds, key=lambda x: -x["shared_pathway_links"])[:200]).to_csv(outdir / "m4_link_prediction.csv", index=False)
    # structural evidence for top targets by Open Targets association
    ranked = (ot.sort_values("open_targets_score", ascending=False)["target"].tolist() if not ot.empty else [])
    ranked = [g for g in ranked if g in set(genes)][:top_structure] or genes[:top_structure]
    srows = []
    for g in ranked:
        acc = uni.get(g) or uniprot_accession(cache, g)
        af, plddt = alphafold_available(cache, acc) if acc else (None, None)
        pdb, pdbid = pdb_available(cache, acc) if acc else (None, None)
        srows.append({"target": g, "uniprot": acc, "alphafold_model": af, "alphafold_mean_plddt": plddt,
                      "experimental_pdb": pdb, "example_pdb_id": pdbid,
                      "docking_status": "pending_requires_vina_gnina_and_ligand_prep"})
    pd.DataFrame(srows).to_csv(outdir / "m4_structure_evidence.csv", index=False)
    pending = pd.DataFrame([
        {"step": "PyKEEN/R-GCN/HGT link prediction with time/book split", "requires": "full heterogeneous KG + trained embedding", "status": "scaffolded_meta_path_only"},
        {"step": "GNINA/AutoDock Vina cross-docking", "requires": "docking binaries + prepared receptors/ligands", "status": "pending_binaries"},
        {"step": "MD + MM/GBSA", "requires": "GROMACS/AMBER + force fields", "status": "pending_binaries"},
        {"step": "PrimeKG/BindingDB integration", "requires": "downloaded PrimeKG CSV + BindingDB TSV", "status": "pending_inputs"},
    ])
    pending.to_csv(outdir / "m4_pending_steps.csv", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--cache", default=".cache/public_api")
    ap.add_argument("--modules", default="1,2,3,4")
    args = ap.parse_args()
    tab = Path(args.outdir) / "tables"; tab.mkdir(parents=True, exist_ok=True)
    cache = ApiCache(Path(args.cache))
    mods = set(args.modules.split(","))
    if "1" in mods: module1_causal_pharmacology(tab, cache); print("Module 1 (genetics-anchored causal pharmacology) done")
    if "2" in mods: module2_niche_localisation(tab); print("Module 2 (single-cell/spatial niche localisation) done")
    if "3" in mods: module3_perturbation_reversal(tab, cache); print("Module 3 (perturbation reversal) done")
    if "4" in mods: module4_knowledge_graph(tab, cache); print("Module 4 (knowledge graph + structure) done")
    print(f"High-order modules written to {tab}")


if __name__ == "__main__":
    main()
