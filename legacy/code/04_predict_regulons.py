"""
Step 4 — Predict regulons via genomic-context / operon analysis

Reproduction code for artifact: predicted_regulons.csv
Artifact version_id: 3bf6113c-1839-4276-9870-f4c5a976504c
Conda environment: python

This is the lineage-tracked code that produced the step's primary output.
Input artifacts are referenced by {{artifact:VERSION_ID}} markers, which the
Claude Science runtime resolves to local file paths; outside that runtime,
replace them with paths to the corresponding CSV/FASTA inputs.
"""

import pandas as pd
import numpy as np
import re
import os
import gzip
import json
import urllib.request
import time
from collections import defaultdict
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

# ── helpers ────────────────────────────────────────────────────────────────────

def iter_fasta(path):
    hdr = None; seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if hdr is not None: yield hdr, "".join(seq)
                hdr = line[1:].split()[0]; seq = []
            else: seq.append(line.strip())
    if hdr is not None: yield hdr, "".join(seq)


def parse_domtbl(fn):
    rows = []
    with open(fn) as f:
        for line in f:
            if line.startswith("#"): continue
            p = line.split()
            if len(p) < 23: continue
            rows.append(dict(
                target=p[0], tacc=p[1], tlen=int(p[2]),
                query=p[3], qacc=p[4], qlen=int(p[5]),
                full_e=float(p[6]), full_score=float(p[7]),
                dom_e=float(p[11]), dom_score=float(p[13]),
                hmm_from=int(p[15]), hmm_to=int(p[16]),
                ali_from=int(p[17]), ali_to=int(p[18]),
            ))
    return pd.DataFrame(rows)


def parse_ft(path):
    rows = []
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if not p or p[0] != "CDS": continue
            try:
                rows.append((p[idx["genomic_accession"]], int(p[idx["start"]]), int(p[idx["end"]]),
                             p[idx["strand"]], p[idx["product_accession"]], p[idx["name"]], p[idx["locus_tag"]]))
            except: continue
    return rows


# ── load project artifacts ─────────────────────────────────────────────────────

# NOTE: host object assumed available in execution environment
A = {a["filename"]: a["latest_version_id"] for a in host.artifacts(project_id="proj_aee37eda121a", limit=400)["artifacts"]}

allg = pd.read_csv(host.artifact_path(A["archaea_all_genomes.csv"]))
per  = pd.read_csv(host.artifact_path(A["regulator_per_genome_full.csv"]))

# ── build discovery panel ──────────────────────────────────────────────────────

rng_panel = np.random.default_rng(42)
per = per.merge(allg[["accession", "ftp_protein_faa", "assembly_level"]], on="accession", how="left")

take_all = ["Korarchaeota (TACK)", "Asgardarchaeota", "DPANN", "Other deep TACK", "Thaumarchaeota (TACK)"]
cap = {"Euryarchaeota": 140, "Crenarchaeota (TACK)": 70}

parts = []
for grp, sub in per.groupby("major_group"):
    if grp in cap:
        n = cap[grp]
        if len(sub) <= n:
            parts.append(sub)
        else:
            frac = n / len(sub)
            s = sub.groupby("subgroup", group_keys=False).apply(
                lambda d: d.sample(max(1, round(len(d) * frac)), random_state=42) if len(d) > 1 else d)
            if len(s) > n: s = s.sample(n, random_state=42)
            parts.append(s)
    else:
        parts.append(sub)

panel = pd.concat(parts).drop_duplicates("accession").reset_index(drop=True)
panel_out = panel[["accession", "species", "genus", "major_group", "subgroup", "n_proteins", "assembly_level", "ftp_protein_faa"]].copy()
panel_out = panel_out.rename(columns={"n_proteins": "proteome_size"})

# ── fetch proteomes ────────────────────────────────────────────────────────────

os.makedirs("proteomes", exist_ok=True)
urls = panel_out[["accession", "ftp_protein_faa"]].values.tolist()

def fix_url(u):
    u = u.replace("//_protein.faa.gz", "/_protein.faa.gz")
    base = u.rsplit("/_protein.faa.gz", 1)[0]
    asm_dir = base.rstrip("/").rsplit("/", 1)[1]
    return f"{base}/{asm_dir}_protein.faa.gz"

fixed = [(acc, fix_url(u)) for acc, u in urls]

ok = 0; fail = []
for i, (acc, u) in enumerate(fixed):
    out = f"proteomes/{acc}.faa.gz"
    if os.path.exists(out) and os.path.getsize(out) > 0:
        ok += 1; continue
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(u), timeout=90) as r:
                data = r.read()
            with open(out, "wb") as f: f.write(data)
            ok += 1; break
        except Exception as e:
            if attempt == 2: fail.append((acc, str(e)))
            else: time.sleep(2)
    if (i + 1) % 50 == 0:
        print(f"{i+1}/{len(fixed)} ok={ok} fail={len(fail)}", flush=True)

# ── concatenate panel proteomes ────────────────────────────────────────────────

os.makedirs("work", exist_ok=True)
with open("work/panel_all.faa", "w") as out_faa:
    for f in sorted(os.listdir("proteomes")):
        if not f.endswith(".faa.gz"): continue
        acc = f[:-7]
        with gzip.open(f"proteomes/{f}", "rt") as gz:
            for line in gz:
                if line.startswith(">"):
                    out_faa.write(f">{acc}|{line[1:]}")
                else:
                    out_faa.write(line)

# ── HMM searches (assumed already run, load results) ──────────────────────────
# These require hmmsearch to have been run externally:
#   hmmsearch --cpu 8 --cut_ga --noali --domtblout work/beyond_known_domtbl.txt beyond_known.hmm work/panel_all.faa
#   hmmsearch --cpu 8 --cut_ga --noali --domtblout work/known_ga_domtbl.txt known19.hmm work/panel_all.faa
#   hmmsearch --cpu 4 -E 1 --noali --domtblout work/known_domtbl.txt known19.hmm work/panel_all.faa

bk = parse_domtbl("work/beyond_known_domtbl.txt")
kn = parse_domtbl("work/known_domtbl.txt")

bk["genome"] = bk["target"].str.split("|").str[0]
bk["protein"] = bk["target"].str.split("|").str[1]
kn["genome"]  = kn["target"].str.split("|").str[0]
kn["protein"] = kn["target"].str.split("|").str[1]

# ── reciprocal exclusion ───────────────────────────────────────────────────────

kn_best   = kn.sort_values("dom_score", ascending=False).groupby("target").first()
kn_score  = kn_best["dom_score"].to_dict()
kn_region = {t: (r.ali_from, r.ali_to) for t, r in kn_best.iterrows()}

def overlap(a, b):
    lo = max(a[0], b[0]); hi = min(a[1], b[1]); return max(0, hi - lo)

def orphan_flag(row):
    t = row["target"]
    if t not in kn_score: return True
    ov = overlap((row.ali_from, row.ali_to), kn_region[t])
    dom_len = row.ali_to - row.ali_from
    frac = ov / max(1, dom_len)
    if frac > 0.5 and kn_score[t] >= row.dom_score - 5:
        return False
    return True

bk["orphan"] = bk.apply(orphan_flag, axis=1)
orph = bk[bk["orphan"]].copy()

def fam_class(q):
    ql = q.lower()
    if ql.startswith("duf") or "duf" in ql: return "uncharacterized_DUF"
    if re.match(r'^(hvo_|mj|pf|saci|ta\d)', ql) or q in ("FilR1_N", "NirdL-like_HTH", "WHD_Lhr", "HVO_A0261_N"):
        return "uncharacterized_cryptic"
    return "characterized_beyond19"

orph["fam_class"] = orph["query"].map(fam_class)

orph_best = orph.sort_values("dom_score", ascending=False).groupby("target").first()

kga = parse_domtbl("work/known_ga_domtbl.txt")
kga["genome"] = kga["target"].str.split("|").str[0]
known_prot = set(kga["target"])

tier = {}
for t, r in orph_best.iterrows():
    if t in known_prot: continue
    tier[t] = "uncharacterized" if r["fam_class"].startswith("uncharacterized") else "characterized_beyond19"
for t in known_prot:
    tier[t] = "known19"

tier_df = pd.DataFrame([{"target": t, "tier": v, "genome": t.split("|")[0]} for t, v in tier.items()])
g2grp = per.set_index("accession")["major_group"].to_dict()
tier_df["group"] = tier_df["genome"].map(g2grp)

# ── clustering (mmseqs2 assumed run externally) ────────────────────────────────
# mmseqs createdb work/orphan_proteins.faa work/mm/orphDB
# mmseqs cluster work/mm/orphDB work/mm/clu work/mm/tmp --min-seq-id 0.30 -c 0.6 --cov-mode 0 --cluster-mode 1 --threads 8
# mmseqs createtsv work/mm/orphDB work/mm/orphDB work/mm/clu work/mm/clusters.tsv

# Extract orphan proteins
orphan_ids = set(orph["target"].unique())
with open("work/orphan_proteins.faa", "w") as out:
    for hdr, seq in iter_fasta("work/panel_all.faa"):
        if hdr in orphan_ids:
            out.write(f">{hdr}\n{seq}\n")

plen = {}
for hdr, seq in iter_fasta("work/orphan_proteins.faa"):
    plen[hdr] = len(seq)

clu = pd.read_csv("work/mm/clusters.tsv", sep="\t", header=None, names=["rep", "member"])
clu["genome"] = clu["member"].str.split("|").str[0]
clu["group"]  = clu["genome"].map(g2grp)
clu["len"]    = clu["member"].map(plen)
mem_fam   = orph_best["query"].to_dict()
mem_class = orph_best["fam_class"].to_dict()
clu["pfam"]       = clu["member"].map(mem_fam)
clu["pfam_class"] = clu["member"].map(mem_class)

def mode_or_mixed(s):
    vc = s.value_counts()
    return vc.index[0] if len(vc) else "NA"

csum = clu.groupby("rep").agg(
    n_members=("member", "size"),
    n_genomes=("genome", "nunique"),
    n_phyla=("group", "nunique"),
    med_len=("len", "median"),
    frac_under250=("len", lambda x: (x < 250).mean()),
    dom_pfam=("pfam", mode_or_mixed),
    dom_class=("pfam_class", mode_or_mixed),
).reset_index()

csum["pass_filter"] = (
    (csum["med_len"] < 250) &
    (csum["frac_under250"] >= 0.6) &
    (csum["n_genomes"] >= 5) &
    (csum["n_phyla"] >= 1)
)
passed = csum[csum["pass_filter"]].sort_values("n_members", ascending=False).reset_index(drop=True)

NON_TF = {
 "eIF-5_eIF-2B", "Ribosomal_S19e", "TFIIE_alpha", "SRP19", "TFIIE_beta", "RNA_pol_Rpb", "S10_plectin",
 "Ribosomal_L29", "Ribosomal_S13", "tRNA_bind", "RuvA_N", "RuvA_C", "DNA_pol", "Rad51", "RecA", "DEAD",
 "Helicase_C", "ResIII", "HNH", "GIY-YIG", "Homeodomain", "Pou", "Ets", "Fork_head", "CBFD_NFYB_HMF",
 "Myb_DNA-binding", "Myb_DNA-bind", "HTH_psq", "zf-", "Integrase", "rve", "Recombinase", "Resolvase",
 "Terminase", "Phage", "DNA_binding_1", "Sigma54_activat", "Sigma70_r2", "Sigma70_r3", "Sigma70_r1",
 "AAA", "Cdc6_C", "Mor", "CsrA", "KH", "OB_", "RRM", "dsrm", "PAZ", "Piwi",
}

def is_non_tf(pf):
    if pd.isna(pf): return False
    for pat in NON_TF:
        if pf == pat or pf.startswith(pat): return True
    return False

passed2 = passed[~passed["dom_pfam"].map(is_non_tf)].copy().reset_index(drop=True)
passed2 = passed2.sort_values("n_members", ascending=False).reset_index(drop=True)
passed2.insert(0, "family_id", [f"nArcR-{i+1}" for i in range(len(passed2))])
passed2["novelty"] = np.where(passed2["dom_class"].str.startswith("uncharacterized"),
                              "uncharacterized", "beyond19_characterized")

rep2fam = dict(zip(passed2["rep"], passed2["family_id"]))
clu_pass = clu[clu["rep"].isin(rep2fam)].copy()
clu_pass["family_id"] = clu_pass["rep"].map(rep2fam)

cat_out = passed2[["family_id", "n_members", "n_genomes", "n_phyla", "med_len", "frac_under250", "dom_pfam", "dom_class", "novelty", "rep"]].copy()
cat_out = cat_out.rename(columns={"dom_pfam": "dominant_pfam", "med_len": "median_length_aa", "rep": "representative_protein"})

fam_genomes = sorted(clu_pass["genome"].unique())

# ── fetch feature tables ───────────────────────────────────────────────────────

panel_urls = dict(zip(panel_out["accession"], panel_out["ftp_protein_faa"]))
ft_urls = []
for acc in fam_genomes:
    u = panel_urls.get(acc)
    if not u: continue
    base = u.rsplit("/", 1)[0]
    asm  = base.rsplit("/", 1)[1]
    ft_urls.append((acc, f"{base}/{asm}_feature_table.txt.gz"))

ft_fixed = []
for acc, u in ft_urls:
    u = u.replace("//_feature_table.txt.gz", "/_feature_table.txt.gz")
    base = u.rsplit("/_feature_table.txt.gz", 1)[0]
    asm = base.rstrip("/").rsplit("/", 1)[1]
    ft_fixed.append((acc, f"{base}/{asm}_feature_table.txt.gz"))

os.makedirs("feature_tables", exist_ok=True)
ok = 0; fail = []
for i, (acc, u) in enumerate(ft_fixed):
    out_path = f"feature_tables/{acc}.ft.txt.gz"
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        ok += 1; continue
    for a in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(u), timeout=90) as r:
                data = r.read()
            open(out_path, "wb").write(data); ok += 1; break
        except Exception as e:
            if a == 2: fail.append((acc, str(e)[:80]))
            else: time.sleep(2)
    if (i + 1) % 80 == 0:
        print(f"{i+1}/{len(ft_fixed)} ok={ok} fail={len(fail)}", flush=True)

# ── build neighborhood table ───────────────────────────────────────────────────

member_by_genome = clu_pass.groupby("genome")["member"].apply(
    lambda s: set(t.split("|")[1] for t in s)).to_dict()
member_fam = dict(zip(
    clu_pass["member"].str.split("|").str[1] + "@" + clu_pass["genome"],
    clu_pass["family_id"]))

neigh_rows = []
missing = 0
genomes_done = 0
for acc in fam_genomes:
    gzp = f"feature_tables/{acc}.ft.txt.gz"
    if not os.path.exists(gzp): missing += 1; continue
    cds = parse_ft(gzp)
    if not cds: continue
    bycontig = defaultdict(list)
    for r in cds: bycontig[r[0]].append(r)
    tgt = member_by_genome[acc]
    for contig, lst in bycontig.items():
        lst.sort(key=lambda x: x[1])
        for i, r in enumerate(lst):
            if r[4] in tgt:
                fam = member_fam.get(r[4] + "@" + acc)
                for off in range(-5, 6):
                    j = i + off
                    if j < 0 or j >= len(lst) or off == 0: continue
                    nb = lst[j]
                    dist = nb[1] - r[2] if off > 0 else r[1] - nb[2]
                    neigh_rows.append((fam, acc, r[4], off, nb[4], nb[5], int(nb[3] == r[3]), dist))
    genomes_done += 1

neigh = pd.DataFrame(neigh_rows, columns=["family_id", "genome", "member_prod", "offset", "nb_prod", "nb_name", "same_strand", "dist"])

# ── functional categorisation ──────────────────────────────────────────────────

CATS = {
 "transporter_ABC":         r"ABC transporter|ATP-binding|permease",
 "transporter_MFS":         r"MFS transporter|major facilitator",
 "transport_other":         r"transporter|symporter|antiporter|channel|efflux",
 "TA_toxin":                r"toxin|VapC|antitoxin|HicA|HicB|RelE|ParE|MazF|PIN ",
 "stress_USP":              r"universal stress",
 "methyltransferase":       r"methyltransferase|SAM-dependent",
 "acetyltransferase":       r"acetyltransferase|GNAT",
 "oxidoreductase":          r"oxidoreductase|dehydrogenase|reductase|FAD|NAD",
 "metallohydrolase":        r"metallo-hydrolase|MBL fold|hydrolase",
 "peptidase":               r"peptidase|protease",
 "kinase":                  r"kinase|phosphatase",
 "signal_HalOD1":           r"HalOD1|output domain|PAS|GAF|signal",
 "chaperone_proteostasis":  r"chaperone|protease|Clp|Hsp|GroEL|DnaK|degradation|proteasome",
 "nucleotide_cyclase":      r"cyclase|diguanylate|c-di-GMP|adenylate",
 "DNA_repair_replication":  r"repair|helicase|nuclease|recombinase|topoisomerase|replicat",
 "metabolism_biosynth":     r"synthase|synthetase|biosynth|aminotransferase|decarboxylase|dehydratase",
 "restriction_modification":r"restriction|modification|methylase",
 "CRISPR_defense":          r"CRISPR|Cas[0-9]|argonaute|defense",
 "regulator_HTH":           r"helix-turn-helix|winged helix|transcriptional regulator|DNA-binding",
 "hypothetical":            r"hypothetical protein|DUF|uncharacterized",
}

def categorize(name):
    if not isinstance(name, str): return []
    hits = [c for c, pat in CATS.items() if re.search(pat, name, re.I)]
    return hits or ["other"]

neigh["cats"] = neigh["nb_name"].map(categorize)
neigh_l = neigh.explode("cats")

# ── build per-genome CDS category pool for null ────────────────────────────────

genome_cds_cats = {}
for acc in fam_genomes:
    gzp = f"feature_tables/{acc}.ft.txt.gz"
    if not os.path.exists(gzp): continue
    cds = parse_ft(gzp)
    cc = []
    for r in cds:
        cc.extend(categorize(r[5]))
    genome_cds_cats[acc] = cc

# ── enrichment tests ───────────────────────────────────────────────────────────

results = []
rng = np.random.default_rng(42)
for fam, sub in neigh_l.groupby("family_id"):
    n = len(sub)
    obs = sub["cats"].value_counts()
    fam_genomes_here = neigh[neigh["family_id"] == fam]["genome"].unique()
    pool = []
    for g in fam_genomes_here:
        pool.extend(genome_cds_cats.get(g, []))
    if not pool: continue
    pool = np.array(pool)
    null_draw = rng.choice(pool, size=min(n * 20, 200000), replace=True)
    null_freq = pd.Series(null_draw).value_counts(normalize=True)
    for cat, k in obs.items():
        if cat in ("other", "hypothetical"): continue
        p0 = null_freq.get(cat, 0.5 / len(null_draw))
        if p0 <= 0: p0 = 0.5 / len(null_draw)
        pv = binomtest(int(k), n, p0, alternative="greater").pvalue
        results.append((fam, cat, int(k), n, k / n, p0, (k / n) / p0 if p0 > 0 else np.nan, pv))

enr = pd.DataFrame(results, columns=["family_id", "category", "obs", "n_slots", "obs_freq", "null_freq", "fold_enrich", "pval"])
enr["padj"] = multipletests(enr["pval"], method="fdr_bh")[1]
enr_sig = enr[(enr["padj"] < 0.05) & (enr["fold_enrich"] >= 1.5) & (enr["obs"] >= 5)].copy()

# ── build predicted_regulons.csv ───────────────────────────────────────────────

cat_meta = cat_out.set_index("family_id")
reg_rows = []
for fam in cat_out["family_id"]:
    sig = enr_sig[enr_sig["family_id"] == fam].sort_values("padj")
    top = sig.head(3)
    preds = "; ".join(f"{r.category}(FE={r.fold_enrich:.1f},padj={r.padj:.1e})" for r in top.itertuples())
    reg_rows.append({
        "family_id": fam,
        "novelty": cat_meta.loc[fam, "novelty"],
        "dominant_pfam": cat_meta.loc[fam, "dominant_pfam"],
        "n_members": cat_meta.loc[fam, "n_members"],
        "n_genomes": cat_meta.loc[fam, "n_genomes"],
        "n_phyla": cat_meta.loc[fam, "n_phyla"],
        "n_sig_categories": len(sig),
        "predicted_regulon_context": preds if preds else "no_significant_enrichment",
        "top_category": top.iloc[0]["category"] if len(top) else "none",
    })

predreg = pd.DataFrame(reg_rows)
predreg.to_csv("predicted_regulons.csv", index=False)
enr_sig.to_csv("neighbor_enrichment_significant.csv", index=False)
print("predicted_regulons.csv written:", len(predreg))