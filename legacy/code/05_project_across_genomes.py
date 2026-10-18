"""
Step 5 — Project novel families across all genomes + MAGs

Reproduction code for artifact: presence_prevalence.csv
Artifact version_id: ff6f2137-0a9b-4a28-92eb-6164e6a1d0c8
Conda environment: python

This is the lineage-tracked code that produced the step's primary output.
Input artifacts are referenced by {{artifact:VERSION_ID}} markers, which the
Claude Science runtime resolves to local file paths; outside that runtime,
replace them with paths to the corresponding CSV/FASTA inputs.
"""

import pandas as pd
import numpy as np
import subprocess
import os

# Load data
pt = pd.read_csv('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/1c109218-69f4-46df-860a-cbb2654c87d1/v32ac8bc7_projection_targets.csv')
g2grp = {}
for _, r in pt.iterrows():
    g2grp.setdefault(str(r["accession"]), r["group"])

panel = pd.read_csv('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/7e6f13ca-e7d8-40e9-87ab-f68e21a06707/v925a23ff_discovery_panel.csv')
for _, r in panel.iterrows():
    g2grp.setdefault(str(r["accession"]), r["major_group"])

def canon(g):
    if not isinstance(g, str): return "Unclassified"
    g = g.strip()
    m = {"Euryarchaeota": "Euryarchaeota", "Crenarchaeota": "Crenarchaeota_TACK", "Crenarchaeota_TACK": "Crenarchaeota_TACK",
         "Thaumarchaeota": "Thaumarchaeota_TACK", "Thaumarchaeota_TACK": "Thaumarchaeota_TACK", "DPANN": "DPANN",
         "Asgardarchaeota": "Asgardarchaeota", "Asgard": "Asgardarchaeota", "Korarchaeota": "Korarchaeota_TACK",
         "Korarchaeota_TACK": "Korarchaeota_TACK", "Other_deep_TACK": "Other_deep_TACK"}
    return m.get(g, g if g in set(m.values()) else "Unclassified")

h = pd.read_csv('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/7e41811b-6038-4df6-949c-50d25b338085/v3e05e2ad_diamond_hits.tsv', sep="\t", header=None,
    names=["qseqid", "sseqid", "pident", "length", "qcov", "scov", "evalue", "bits"])
h["genome"] = h["qseqid"].str.split("|").str[0]
h["family"] = h["sseqid"].str.split("__").str[0]
h["phylum"] = h["genome"].map(lambda g: canon(g2grp.get(g)))

phyla_order = ["Euryarchaeota", "Crenarchaeota_TACK", "Thaumarchaeota_TACK", "Other_deep_TACK", "DPANN", "Asgardarchaeota", "Korarchaeota_TACK"]
fams = sorted(h["family"].unique(), key=lambda x: int(x.split("-")[1]))

# Build scanned_ph from the hit genomes (approximation without the full short.faa)
scanned_ph = h.drop_duplicates("genome").set_index("genome")["phylum"]

cn = h.groupby(["family", "genome"]).size().unstack(fill_value=0).reindex(index=fams, fill_value=0)
pa = (cn > 0).astype(int)

prev = pd.DataFrame(index=fams, columns=phyla_order, dtype=float)
gph = scanned_ph
for ph in phyla_order:
    ph_genomes = set(gph[gph == ph].index)
    d = len(ph_genomes)
    if d == 0:
        prev[ph] = 0.0
        continue
    present_cols = [g for g in pa.columns if g in ph_genomes]
    sub = pa[present_cols] if present_cols else pd.DataFrame(index=pa.index)
    prev[ph] = (sub.sum(axis=1) / d) if present_cols else 0.0
prev = prev.fillna(0.0)

prev_out = prev.copy()
prev_out.insert(0, "family", prev_out.index)
prev_out.to_csv("presence_prevalence.csv", index=False)