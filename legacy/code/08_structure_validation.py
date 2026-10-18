"""
Step 8 — Structure-based validation on a representative shortlist

Reproduction code for artifact: structure_validation.csv
Artifact version_id: bc7689ce-031f-420c-8128-e5cc2eeacfb4
Conda environment: python

This is the lineage-tracked code that produced the step's primary output.
Input artifacts are referenced by {{artifact:VERSION_ID}} markers, which the
Claude Science runtime resolves to local file paths; outside that runtime,
replace them with paths to the corresponding CSV/FASTA inputs.
"""

import pandas as pd
import numpy as np
import json
import os
import glob
import requests
import time

# Load catalog and membership
cat = pd.read_csv('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/6ca67a0c-f30f-4889-8af2-574a9937463c/v20b91b91_novel_family_catalog.csv')
mem = pd.read_csv('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/2062e495-9137-4c58-90fa-ddb83c8341c2/v0c4ca4a6_novel_family_membership.csv')

# Load AF mapping
mapping = json.load(open('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/b7118a5b-f8a4-4b3c-9dae-4ad27148a7a9/v6d34951a_af_mapping.json'))

# Fetch AlphaFold DB PDB models
os.makedirs("val/af_pdb", exist_ok=True)
sess = requests.Session()
ok = 0
fail = []
t0 = time.time()
for i, (fam, (acc, up, afdb)) in enumerate(mapping.items()):
    fp = f"val/af_pdb/{fam}.pdb"
    if os.path.exists(fp) and os.path.getsize(fp) > 1000:
        ok += 1
        continue
    try:
        a = sess.get(f"https://alphafold.ebi.ac.uk/api/prediction/{afdb}", timeout=20)
        if a.status_code == 200 and a.json():
            url = a.json()[0]["pdbUrl"]
            r = sess.get(url, timeout=30)
            if r.status_code == 200 and len(r.text) > 1000:
                open(fp, "w").write(r.text)
                ok += 1
                continue
    except Exception:
        pass
    fail.append(fam)
    if (i + 1) % 50 == 0:
        print(f"{i+1}/{len(mapping)} ok={ok} {round(time.time()-t0)}s")
print("DONE ok=", ok, "fail=", len(fail))

import biotite.structure.io.pdb as pdb
import biotite.structure as struc

def analyze(fp):
    f = pdb.PDBFile.read(fp)
    arr = f.get_structure(model=1, extra_fields=["b_factor"])
    arr = arr[struc.filter_amino_acids(arr)]
    ca = arr[arr.atom_name == "CA"]
    plddt = float(np.mean(ca.b_factor)) if len(ca) else np.nan
    nres = len(ca)
    try:
        sse = struc.annotate_sse(arr)
        h = float(np.mean(sse == "a"))
        s = float(np.mean(sse == "b"))
        c = float(np.mean(sse == "c"))
    except Exception:
        h = s = c = np.nan
    return nres, plddt, h, s, c

rows = []
for fp in sorted(glob.glob("val/af_pdb/*.pdb")):
    fam = os.path.basename(fp)[:-4]
    try:
        nres, plddt, h, s, c = analyze(fp)
        rows.append((fam, nres, round(plddt, 1), round(h, 3), round(s, 3), round(c, 3)))
    except Exception as e:
        rows.append((fam, None, None, None, None, None))
S = pd.DataFrame(rows, columns=["family_id", "n_res", "mean_plddt", "frac_helix", "frac_sheet", "frac_coil"])

S = S.merge(cat[["family_id", "novelty", "dominant_pfam", "dom_class", "n_phyla"]], on="family_id", how="left")
S["uniprot"] = S.family_id.map(lambda f: mapping[f][1] if f in mapping else None)

def verdict(r):
    conf = r.mean_plddt >= 70
    helix_rich = r.frac_helix >= 0.30
    helix_dom = r.frac_helix >= r.frac_sheet
    if conf and helix_rich and helix_dom:
        return "PASS: confident DNA-binding-like fold"
    if r.mean_plddt >= 50 and helix_dom and r.frac_helix >= 0.20:
        return "PASS: helix-dominant (moderate conf)"
    return "REVIEW"

S["structure_verdict"] = S.apply(verdict, axis=1)
S.to_csv("structure_validation.csv", index=False)