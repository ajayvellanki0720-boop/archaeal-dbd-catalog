"""
Step 7 — Domain-level validation of family representatives

Reproduction code for artifact: domain_validation.csv
Artifact version_id: f88d8633-33a2-48e1-9e6b-5d1d81fc64bb
Conda environment: python

This is the lineage-tracked code that produced the step's primary output.
Input artifacts are referenced by {{artifact:VERSION_ID}} markers, which the
Claude Science runtime resolves to local file paths; outside that runtime,
replace them with paths to the corresponding CSV/FASTA inputs.
"""

import pandas as pd

cat = pd.read_csv('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/6ca67a0c-f30f-4889-8af2-574a9937463c/v20b91b91_novel_family_catalog.csv')
mem = pd.read_csv('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/2062e495-9137-4c58-90fa-ddb83c8341c2/v0c4ca4a6_novel_family_membership.csv')
sc = pd.read_csv('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/4e19d9f3-6caf-46d3-af52-725aea6b3daf/v93ea0b71_self_consistency.csv')

k19 = {"HTH_1","HTH_AraC","Crp","GntR","TetR_N","HTH_5","AsnC_trans_reg","MarR","HTH_3","RHH_1","HTH_6","TrmB","LysR_substrate","HTH_10","Sigma70_r4_2","HTH_17","Acetyltransf_3","HTH_23","MerR_1","HTH_31"}

selfbest = sc.groupby("true_fam")["self_best"].mean()
arch = mem.groupby("family_id").apply(lambda g: g["pfam"].fillna("none").value_counts().iloc[0]/len(g), include_groups=False)

v = cat[["family_id","dominant_pfam","dom_class","novelty","n_members","n_genomes","n_phyla"]].copy()
v["self_best_frac"] = v.family_id.map(selfbest)
v["arch_homogeneity"] = v.family_id.map(arch)
v["dom_pfam_is_known19"] = v.dominant_pfam.isin(k19)
v["domain_validation"] = ["PASS: coherent + homogeneous" if(a>=0.99 and b>=0.8) else("PASS: coherent" if a>=0.5 else "REVIEW: mixed") for a,b in zip(v.self_best_frac, v.arch_homogeneity)]

R = pd.read_csv('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/544107ea-7539-431e-9491-1db65f75deb6/v59a91798_reps_vs_pfamA.domtbl', sep='\t', header=None,
                names=["rep","pfam","seqE","seqscore"])
R["fam"] = R["rep"].str.split("__").str[0]
best = R.sort_values("seqscore", ascending=False).groupby("rep").first().reset_index()
best["fam"] = best["rep"].str.split("__").str[0]

fam_any = R.groupby("fam")["rep"].nunique()
fam_bestpfam = best.groupby("fam")["pfam"].apply(lambda s: s.mode().iloc[0] if len(s) else None)

cat2 = cat.set_index("family_id")
res = []
for fam in cat.family_id:
    nrep_hit = int(fam_any.get(fam, 0))
    bp = fam_bestpfam.get(fam, None)
    nov = cat2.loc[fam, "novelty"]
    dompf = cat2.loc[fam, "dominant_pfam"]
    is_duf = bool(bp and ("DUF" in str(bp) or "UPF" in str(bp)))
    res.append((fam, nov, dompf, nrep_hit, bp, is_duf))
Rdf = pd.DataFrame(res, columns=["family_id","novelty","catalog_pfam","n_reps_with_pfam_hit","best_pfam_independent","best_is_DUF"])

v = v.merge(Rdf[["family_id","n_reps_with_pfam_hit","best_pfam_independent","best_is_DUF"]], on="family_id", how="left")
v["pfam_recovers_catalog"] = v["best_pfam_independent"] == v["dominant_pfam"]

v.to_csv("domain_validation.csv", index=False)