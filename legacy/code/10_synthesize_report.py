"""
Step 10 — Synthesize map, predictions, evolution, and validation

Reproduction code for artifact: MANUSCRIPT.md
Artifact version_id: 87ffb636-2482-4722-b11c-ffd9c7c6af3e
Conda environment: python

This is the lineage-tracked code that produced the step's primary output.
Input artifacts are referenced by {{artifact:VERSION_ID}} markers, which the
Claude Science runtime resolves to local file paths; outside that runtime,
replace them with paths to the corresponding CSV/FASTA inputs.
"""

# [lineage] reconstruction failed validation (max_tokens) — raw producing cell below; see Execution Log for full trace
p="MANUSCRIPT.md"
s=open(p).read()
# add software paragraph to 4.8 reproducibility
old="Random seeds are fixed (42) throughout."
new=("Random seeds are fixed (42) throughout. Data handling, statistics, and visualization used "
     "NumPy [32], SciPy [33], pandas [34], statsmodels [35], Matplotlib [36], and Biopython [37]; "
     "differential expression used edgeR-style CPM normalization and dispersion estimation.")
assert old in s
s=s.replace(old,new,1)

refs=open("val/references_section.md").read()
# append references + statements before nothing (end of doc)
tail=f"""

---

## References

{refs}

---

## Declarations

**Data and code availability.** All datasets analyzed here are public: archaeal genome assemblies from NCBI GenBank/RefSeq, Pfam-A profile HMMs, AlphaFold Database structure models, and *Haloferax volcanii* RNA-seq series GSE293312 (NCBI GEO). Every analytical stage is backed by a lineage-tracked project artifact (Section 5); random seeds are fixed at 42 throughout for reproducibility.

**Competing interests.** The authors declare no competing interests.

**Funding.** Not applicable.

**Author contributions.** Conceptualization, analysis, and writing were performed within the Claude Science computational environment (see Section 6).
"""
s=s.rstrip()+tail
open(p,"w").write(s)
print("bytes:",len(s))
print(s[-500:])
