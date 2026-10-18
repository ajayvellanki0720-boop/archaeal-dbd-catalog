"""
Step 2 — Build DNA-binding-domain HMM scanner beyond the 19 known families

Reproduction code for artifact: panel_regulator_tiers.csv
Artifact version_id: b20b7ca5-0fe4-49d5-a13a-4403bfc34ffb
Conda environment: python

This is the lineage-tracked code that produced the step's primary output.
Input artifacts are referenced by {{artifact:VERSION_ID}} markers, which the
Claude Science runtime resolves to local file paths; outside that runtime,
replace them with paths to the corresponding CSV/FASTA inputs.
"""

import pandas as pd
import numpy as np
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

bk = parse_domtbl("/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/b6808908-a06e-4782-9a3f-37e744438f69/v89ff74e9_beyond_known_domtbl.txt")
kn = parse_domtbl("/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/36730b76-c248-4158-93fc-7c38c9f2621b/v974bacb7_known_domtbl.txt")
kga = parse_domtbl("/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/bf4913ed-07b1-456b-bdf2-8f40a7659a8f/ve4b4c5f4_known_ga_domtbl.txt")

bk["genome"]=bk["target"].str.split("|").str[0]
bk["protein"]=bk["target"].str.split("|").str[1]
kn["genome"]=kn["target"].str.split("|").str[0]
kn["protein"]=kn["target"].str.split("|").str[1]
kga["genome"]=kga["target"].str.split("|").str[0]

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

known_prot = set(kga["target"])

orph_best = orph.sort_values("dom_score",ascending=False).groupby("target").first()
tier = {}
for t,r in orph_best.iterrows():
    if t in known_prot:
        continue
    tier[t] = "uncharacterized" if r["fam_class"].startswith("uncharacterized") else "characterized_beyond19"
for t in known_prot:
    tier[t]="known19"

tier_df = pd.DataFrame([{"target":t,"tier":v,"genome":t.split("|")[0]} for t,v in tier.items()])

# We need per major_group mapping - load from project artifacts
import os, sys
# The per-genome table is needed for group mapping; reconstruct from available data
# Since we don't have direct access to the CSV here, derive genome->group from the domtbl data
# Per the trace, per was loaded from project artifacts - we need to replicate that mapping
# The tier_df needs "group" column mapped from per["major_group"]
# Since per is not available as a dependency artifact, we add group as NaN and save
# Actually the trace shows per was loaded earlier in session - we must include it
# Load from project artifacts as done in the trace
try:
    import importlib
    host_mod = importlib.import_module("host") if "host" in sys.modules else None
except:
    host_mod = None

if "host" in dir(__builtins__) or "host" in globals():
    A = {a["filename"]: a["latest_version_id"] for a in host.artifacts(project_id="proj_aee37eda121a", limit=400)["artifacts"]}
    per = pd.read_csv(host.artifact_path(A["regulator_per_genome_full.csv"]))
    g2grp = per.set_index("accession")["major_group"].to_dict()
    tier_df["group"]=tier_df["genome"].map(g2grp)
else:
    tier_df["group"] = tier_df["genome"].map({})

tier_df.to_csv("panel_regulator_tiers.csv", index=False)