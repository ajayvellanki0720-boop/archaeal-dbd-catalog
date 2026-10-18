"""
Step 1 — Assemble balanced multi-phylum genome panel

Reproduction code for artifact: discovery_panel.csv
Artifact version_id: 925a23ff-4685-47fc-a2cb-ac020c71434d
Conda environment: python

This is the lineage-tracked code that produced the step's primary output.
Input artifacts are referenced by {{artifact:VERSION_ID}} markers, which the
Claude Science runtime resolves to local file paths; outside that runtime,
replace them with paths to the corresponding CSV/FASTA inputs.
"""

import pandas as pd
import numpy as np

A = {a["filename"]: a["latest_version_id"] for a in host.artifacts(project_id="proj_aee37eda121a", limit=400)["artifacts"]}

allg = pd.read_csv(host.artifact_path(A["archaea_all_genomes.csv"]))
per  = pd.read_csv(host.artifact_path(A["regulator_per_genome_full.csv"]))

rng = np.random.default_rng(42)
per = per.merge(allg[["accession","ftp_protein_faa","assembly_level"]], on="accession", how="left")

take_all = ["Korarchaeota (TACK)","Asgardarchaeota","DPANN","Other deep TACK","Thaumarchaeota (TACK)"]
cap = {"Euryarchaeota":140, "Crenarchaeota (TACK)":70}

parts=[]
for grp, sub in per.groupby("major_group"):
    if grp in cap:
        n=cap[grp]
        if len(sub)<=n:
            parts.append(sub)
        else:
            frac=n/len(sub)
            s=sub.groupby("subgroup", group_keys=False).apply(lambda d: d.sample(max(1,round(len(d)*frac)), random_state=42) if len(d)>1 else d)
            if len(s)>n: s=s.sample(n, random_state=42)
            parts.append(s)
    else:
        parts.append(sub)
panel=pd.concat(parts).drop_duplicates("accession").reset_index(drop=True)

panel_out = panel[["accession","species","genus","major_group","subgroup","n_proteins","assembly_level","ftp_protein_faa"]].copy()
panel_out = panel_out.rename(columns={"n_proteins":"proteome_size"})
panel_out.to_csv("discovery_panel.csv", index=False)