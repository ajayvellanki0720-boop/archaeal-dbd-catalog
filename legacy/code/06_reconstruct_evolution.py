"""
Step 6 — Reconstruct evolution of top novel families across phyla

Reproduction code for artifact: family_evolution_congruence.csv
Artifact version_id: 1db3f1aa-94df-474e-b1b7-029e67a2e800
Conda environment: python

This is the lineage-tracked code that produced the step's primary output.
Input artifacts are referenced by {{artifact:VERSION_ID}} markers, which the
Claude Science runtime resolves to local file paths; outside that runtime,
replace them with paths to the corresponding CSV/FASTA inputs.
"""

import pandas as pd
import numpy as np
import random
from Bio import Phylo

random.seed(42)

tree_families = ['nArcR-1', 'nArcR-32', 'nArcR-8', 'nArcR-7', 'nArcR-51', 'nArcR-11', 'nArcR-33', 'nArcR-10', 'nArcR-16', 'nArcR-17']

nwk_paths = {
    'nArcR-1': '/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/f360c138-092f-4705-93d4-dd834b7bcbc3/va5a546a9_nArcR-1.nwk',
    'nArcR-10': '/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/28398f7f-559f-449a-a965-3ab37af7a122/vfae3180f_nArcR-10.nwk',
    'nArcR-11': '/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/e90269cc-675b-44f5-b529-6c6aa1fe4ad8/v3d2c2f5a_nArcR-11.nwk',
    'nArcR-16': '/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/6ead2e87-1467-4d19-a1de-830321b69532/va255d0b3_nArcR-16.nwk',
    'nArcR-17': '/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/b4e08e80-9db8-40ae-93a5-d286e347b4f5/v46aa4095_nArcR-17.nwk',
    'nArcR-32': '/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/dfb6f8cf-c9fa-4e27-a70e-5bb907a7fab8/v21de4a89_nArcR-32.nwk',
    'nArcR-33': '/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/92994406-0935-4558-ae21-183fb85057e1/v49342ff9_nArcR-33.nwk',
    'nArcR-51': '/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/d51f5670-3b82-471a-9dba-195055438470/v1380833a_nArcR-51.nwk',
    'nArcR-7': '/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/004e56cd-deda-43b8-8bfa-3356006eb5c0/v48104118_nArcR-7.nwk',
    'nArcR-8': '/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/d469c820-6f12-496d-96fc-6fa1dfc8e765/v19e9935f_nArcR-8.nwk',
}

def congruence(nwk):
    t = Phylo.read(nwk, "newick")
    tips = [l.name for l in t.get_terminals()]
    phs = [name.split("__")[-1] for name in tips]
    same = sum(1 for i in range(len(phs)-1) if phs[i]==phs[i+1])/(len(phs)-1)
    nulls = []
    for _ in range(200):
        sh = phs[:]
        random.shuffle(sh)
        nulls.append(sum(1 for i in range(len(sh)-1) if sh[i]==sh[i+1])/(len(sh)-1))
    nulls = np.array(nulls)
    z = (same - nulls.mean())/(nulls.std()+1e-9)
    p = (np.sum(nulls>=same)+1)/(len(nulls)+1)
    return dict(n_tips=len(tips), n_phyla=len(set(phs)), adj_same=round(same,3),
                null_mean=round(nulls.mean(),3), z=round(z,2), p_monophyly=round(p,4))

rows = []
for fam in tree_families:
    d = congruence(nwk_paths[fam])
    d["family"] = fam
    rows.append(d)

cong = pd.DataFrame(rows)[["family","n_tips","n_phyla","adj_same","null_mean","z","p_monophyly"]]

def interp(r):
    if r["p_monophyly"] > 0.05: return "no phylum structure (single-phylum or panmictic)"
    if r["z"] >= 8: return "strong vertical inheritance (within-phylum clades)"
    if r["z"] >= 4: return "moderate vertical signal with some mixing"
    return "weak structure: ancient/pan-archaeal origin or lateral spread"

cong_out = cong.copy()
cong_out["evolutionary_pattern"] = cong_out.apply(interp, axis=1)

cat = pd.read_csv('/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/artifacts/proj_a956442fed5c/6ca67a0c-f30f-4889-8af2-574a9937463c/v20b91b91_novel_family_catalog.csv')
cong_out = cong_out.merge(cat[["family_id","dominant_pfam","novelty"]], left_on="family", right_on="family_id", how="left").drop(columns="family_id")

cong_out.to_csv("family_evolution_congruence.csv", index=False)