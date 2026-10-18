"""
Step 9 — Transcriptomic validation in model archaea

Reproduction code for artifact: oxidative_response.csv
Artifact version_id: 8a4cd441-2147-4811-ac30-9bffde0cbff5
Conda environment: python

This is the lineage-tracked code that produced the step's primary output.
Input artifacts are referenced by {{artifact:VERSION_ID}} markers, which the
Claude Science runtime resolves to local file paths; outside that runtime,
replace them with paths to the corresponding CSV/FASTA inputs.
"""

import pandas as pd, numpy as np, glob, re, gzip, os
from scipy import stats
from statsmodels.stats.multitest import multipletests

os.makedirs("val/rnaseq", exist_ok=True)

# Load count matrix files
def load_counts(fp):
    df = pd.read_csv(fp, sep="\t", comment="#")
    val = df.columns[-1]
    return df.set_index("Geneid")[val]

files = sorted(glob.glob("val/rnaseq/GSM*_readcount.txt.gz"))
cols = {}
for fp in files:
    m = re.search(r"GSM\d+_(.+?)_RNA_readcount", fp)
    cols[m.group(1)] = load_counts(fp)
M = pd.DataFrame(cols)

M = M[~M.index.duplicated(keep="first")]

# Load HVO novel loci
best = pd.read_csv("val/rnaseq/hvo_novel_loci.csv")
novel_loci = set(best.locus.dropna())

cpm = M / M.sum() * 1e6
wt = ["WTA_cont", "WTB_cont", "WTC_cont"]
oxi = ["WTA_oxi", "WTB_oxi", "WTC_oxi"]

res = []
for g in M.index:
    a = cpm.loc[g, oxi].values.astype(float)
    b = cpm.loc[g, wt].values.astype(float)
    lfc = np.log2((a.mean() + 1) / (b.mean() + 1))
    try:
        t, p = stats.ttest_ind(np.log2(a + 1), np.log2(b + 1))
    except Exception:
        p = np.nan
    res.append((g, lfc, p))

D = pd.DataFrame(res, columns=["locus", "log2FC_oxi_vs_ctrl", "pval"]).set_index("locus")
D["novel"] = D.index.isin(novel_loci)

ok = D.pval.notna().values
padj = np.full(len(D), np.nan)
padj[ok] = multipletests(D.pval.values[ok], method="fdr_bh")[1]
D["padj"] = padj

D.to_csv("val/rnaseq/oxidative_response.csv")