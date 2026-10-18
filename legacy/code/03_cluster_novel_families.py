"""
Step 3 — Cluster orphan candidates into novel regulator families

Reproduction code for artifact: novel_family_catalog.csv
Artifact version_id: 20b91b91-11a2-4dcf-bd1b-8e48daca70d2
Conda environment: archaea-dormancy

This is the lineage-tracked code that produced the step's primary output.
Input artifacts are referenced by {{artifact:VERSION_ID}} markers, which the
Claude Science runtime resolves to local file paths; outside that runtime,
replace them with paths to the corresponding CSV/FASTA inputs.
"""

import pandas as pd
import numpy as np

# Load clusters
clu = pd.read_csv("/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/d83f1ab8-255a-4390-ad3a-2199596aef35/vcd7046e4_clusters.tsv", sep="\t", header=None, names=["rep","member"])

# protein length map
def iter_fasta(path):
    hdr=None; seq=[]
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if hdr is not None: yield hdr, "".join(seq)
                hdr=line[1:].split()[0]; seq=[]
            else: seq.append(line.strip())
    if hdr is not None: yield hdr,"".join(seq)

plen = {}
for hdr,seq in iter_fasta("/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/4fead882-50b4-4e3c-be0e-a8ff7bb2dc2e/va7000a14_orphan_proteins.faa"):
    plen[hdr]=len(seq)

clu["genome"]=clu["member"].str.split("|").str[0]
clu["len"]=clu["member"].map(plen)

# We need group (phylum) mapping - reconstruct from orphan hits context
# The orphan hits and per-genome data need to be reconstructed

import re

def parse_domtbl(fn):
    rows=[]
    with open(fn) as f:
        for line in f:
            if line.startswith("#"): continue
            p=line.split()
            if len(p)<23: continue
            rows.append(dict(
                target=p[0], tacc=p[1], tlen=int(p[2]),
                query=p[3], qacc=p[4], qlen=int(p[5]),
                full_e=float(p[6]), full_score=float(p[7]),
                dom_e=float(p[11]), dom_score=float(p[13]),
                hmm_from=int(p[15]), hmm_to=int(p[16]),
                ali_from=int(p[17]), ali_to=int(p[18]),
            ))
    return pd.DataFrame(rows)

import host

A = {a["filename"]: a["latest_version_id"] for a in host.artifacts(project_id="proj_aee37eda121a", limit=400)["artifacts"]}

per  = pd.read_csv(host.artifact_path(A["regulator_per_genome_full.csv"]))

g2grp = per.set_index("accession")["major_group"].to_dict()
clu["group"]=clu["genome"].map(g2grp)

bk = parse_domtbl("work/beyond_known_domtbl.txt")
kn = parse_domtbl("work/known_domtbl.txt")

bk["genome"]=bk["target"].str.split("|").str[0]
bk["protein"]=bk["target"].str.split("|").str[1]
kn["genome"]=kn["target"].str.split("|").str[0]
kn["protein"]=kn["target"].str.split("|").str[1]

kn_best = kn.sort_values("dom_score", ascending=False).groupby("target").first()
kn_score = kn_best["dom_score"].to_dict()
kn_region = {t:(r.ali_from, r.ali_to) for t,r in kn_best.iterrows()}

def overlap(a,b):
    lo=max(a[0],b[0]); hi=min(a[1],b[1]); return max(0,hi-lo)

def orphan_flag(row):
    t=row["target"]
    if t not in kn_score: return True
    ov = overlap((row.ali_from,row.ali_to), kn_region[t])
    dom_len = row.ali_to-row.ali_from
    frac = ov/max(1,dom_len)
    if frac>0.5 and kn_score[t] >= row.dom_score - 5:
        return False
    return True

bk["orphan"]=bk.apply(orphan_flag, axis=1)
orph = bk[bk["orphan"]].copy()

def fam_class(q):
    ql=q.lower()
    if ql.startswith("duf") or "duf" in ql: return "uncharacterized_DUF"
    if re.match(r'^(hvo_|mj|pf|saci|ta\d)', ql) or q in ("FilR1_N","NirdL-like_HTH","WHD_Lhr","HVO_A0261_N"):
        return "uncharacterized_cryptic"
    return "characterized_beyond19"

orph["fam_class"]=orph["query"].map(fam_class)

orph_best = orph.sort_values("dom_score",ascending=False).groupby("target").first()

mem_fam = orph_best["query"].to_dict()
mem_class = orph_best["fam_class"].to_dict()
clu["pfam"]=clu["member"].map(mem_fam)
clu["pfam_class"]=clu["member"].map(mem_class)

def mode_or_mixed(s):
    vc=s.value_counts()
    return vc.index[0] if len(vc) else "NA"

csum = clu.groupby("rep").agg(
    n_members=("member","size"),
    n_genomes=("genome","nunique"),
    n_phyla=("group","nunique"),
    med_len=("len","median"),
    frac_under250=("len", lambda x:(x<250).mean()),
    dom_pfam=("pfam", mode_or_mixed),
    dom_class=("pfam_class", mode_or_mixed),
).reset_index()

csum["pass_filter"] = (
    (csum["med_len"]<250) &
    (csum["frac_under250"]>=0.6) &
    (csum["n_genomes"]>=5) &
    (csum["n_phyla"]>=1)
)
passed = csum[csum["pass_filter"]].sort_values("n_members",ascending=False).reset_index(drop=True)

NON_TF = {
 "eIF-5_eIF-2B","Ribosomal_S19e","TFIIE_alpha","SRP19","TFIIE_beta","RNA_pol_Rpb","S10_plectin",
 "Ribosomal_L29","Ribosomal_S13","tRNA_bind","RuvA_N","RuvA_C","DNA_pol","Rad51","RecA","DEAD",
 "Helicase_C","ResIII","HNH","GIY-YIG","Homeodomain","Pou","Ets","Fork_head","CBFD_NFYB_HMF",
 "Myb_DNA-binding","Myb_DNA-bind","HTH_psq","zf-","Integrase","rve","Recombinase","Resolvase",
 "Terminase","Phage","DNA_binding_1","Sigma54_activat","Sigma70_r2","Sigma70_r3","Sigma70_r1",
 "AAA","Cdc6_C","Mor","CsrA","KH","OB_","RRM","dsrm","PAZ","Piwi",
}
def is_non_tf(pf):
    if pd.isna(pf): return False
    for pat in NON_TF:
        if pf==pat or pf.startswith(pat): return True
    return False

passed2 = passed[~passed["dom_pfam"].map(is_non_tf)].copy().reset_index(drop=True)

passed2 = passed2.sort_values("n_members", ascending=False).reset_index(drop=True)
passed2.insert(0, "family_id", [f"nArcR-{i+1}" for i in range(len(passed2))])
passed2["novelty"] = np.where(passed2["dom_class"].str.startswith("uncharacterized"),
                              "uncharacterized", "beyond19_characterized")

cat_out = passed2[["family_id","n_members","n_genomes","n_phyla","med_len","frac_under250","dom_pfam","dom_class","novelty","rep"]].copy()
cat_out = cat_out.rename(columns={"dom_pfam":"dominant_pfam","med_len":"median_length_aa","rep":"representative_protein"})
cat_out.to_csv("novel_family_catalog.csv", index=False)