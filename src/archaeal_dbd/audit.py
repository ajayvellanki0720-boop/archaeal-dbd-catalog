from __future__ import annotations

import json
import tarfile
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd


MISSING_CITED_TABLES = (
    "domain_validation.csv",
    "enrichment_block_permutation.csv",
    "enrichment_dependence_corrected.csv",
    "motif_reanalysis_deep.csv",
    "novel_family_distribution.csv",
    "pfam_independent_validation.csv",
    "regulon_motifs.csv",
    "self_consistency.csv",
    "structural_all_hits_expanded.csv",
)


@dataclass(frozen=True)
class Check:
    name: str
    status: str
    observed: object
    expected: object
    detail: str


def _check(name: str, observed: object, expected: object, detail: str = "") -> Check:
    return Check(name, "pass" if observed == expected else "fail", observed, expected, detail)


def run_release_audit(release_dir: Path, output_path: Path | None = None) -> dict:
    data = release_dir / "data"
    tables = release_dir / "tables"
    figures = release_dir / "figures"
    checks: list[Check] = []

    catalog = pd.read_csv(tables / "novel_family_catalog.csv")
    membership = pd.read_csv(data / "novel_family_membership.csv")
    tiers = pd.read_csv(tables / "panel_regulator_tiers.csv")
    accounting = pd.read_csv(tables / "genome_accounting.csv")
    projection_metadata = pd.read_csv(data / "projection_genome_metadata.csv")
    presence = pd.read_csv(data / "presence_absence_matrix.csv", nrows=1)
    copy_number = pd.read_csv(data / "copynumber_matrix.csv", nrows=1)
    motifs = pd.read_csv(tables / "motif_null_analysis.csv")
    structures = pd.read_csv(tables / "structural_novelty_expanded.csv")
    unresolved = pd.read_csv(tables / "uncharacterized_db_comparison.csv")

    checks.extend(
        [
            _check("catalog families", catalog["family_id"].nunique(), 250),
            _check("catalog members", int(catalog["n_members"].sum()), 13_237),
            _check("membership rows", len(membership), 13_237),
            _check("target-clan proteins", len(tiers), 38_153),
            _check("genome accounting rows", len(accounting), 3_526),
            _check("legacy projection metadata rows", len(projection_metadata), 3_192,
                   "This is a legacy intermediate, not the final denominator."),
            _check("legacy presence matrix genomes", len(presence.columns) - 1, 3_192,
                   "A corrected 3,526-genome matrix is not present."),
            _check("legacy copy matrix genomes", len(copy_number.columns) - 1, 3_192,
                   "A corrected 3,526-genome matrix is not present."),
            _check("motifs passing all nulls", int(motifs["passes_all_nulls"].sum()), 11),
            _check("scored structures", len(structures), 184),
            _check("structures matching curated DBD folds", int(structures["matches_TF_055"].sum()), 175),
            _check("unresolved sequence and structure families",
                   int((unresolved["resolution"] == "genuinely_unresolved").sum()), 9),
        ]
    )

    fate_counts = accounting["fate"].value_counts().to_dict()
    expected_fates = {"analyzed_hit_diamond": 3144, "analyzed_hit": 379, "analyzed_zero": 3}
    checks.append(_check("projection accounting fates", fate_counts, expected_fates))

    archive = data / "iqtree_phylogenetics.tar.gz"
    with tarfile.open(archive, "r:gz") as handle:
        real_files = [
            member
            for member in handle.getmembers()
            if member.isfile() and not Path(member.name).name.startswith("._")
        ]
    checks.append(_check("IQ-TREE archive files", len(real_files), 30))

    absent = [name for name in MISSING_CITED_TABLES if not (tables / name).exists()]
    checks.append(_check("manuscript-cited tables absent from release", absent, list(MISSING_CITED_TABLES),
                         "Pass means the audit reproduced the known omission list, not that the release is complete."))

    result = {
        "release_dir": str(release_dir),
        "all_numeric_invariants_pass": all(c.status == "pass" for c in checks),
        "checks": [asdict(c) for c in checks],
        "interpretation": (
            "The deposited summaries are internally checkable, but the 3,526-genome raw matrices, "
            "several cited source tables, and raw reciprocal-exclusion hits are not present."
        ),
    }
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result
