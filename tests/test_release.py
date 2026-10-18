from pathlib import Path

import pandas as pd

from archaeal_dbd import __version__
from archaeal_dbd.audit import run_release_audit
from archaeal_dbd.motif_holdout import auc, dinucleotide_shuffle
from archaeal_dbd.panel import run_panel_sensitivity
from archaeal_dbd.projection import run_independent_projection_resampling
from archaeal_dbd.quality_sensitivity import run_quality_sensitivity
from archaeal_dbd.thresholds import reciprocal_exclusion_grid


ROOT = Path(__file__).resolve().parents[1]


def test_package_version_matches_release():
    assert __version__ == "0.4.0"


def test_release_invariants(tmp_path):
    report = run_release_audit(ROOT / "release", tmp_path / "audit.json")
    assert report["all_numeric_invariants_pass"]


def test_panel_sensitivity_is_deterministic(tmp_path):
    out = run_panel_sensitivity(
        ROOT / "release/tables/discovery_panel.csv",
        ROOT / "release/data/novel_family_membership.csv",
        tmp_path,
        replicates=10,
        seed=10_000,
    )
    summary = pd.read_csv(out.summary).set_index("scheme")
    assert summary.loc["proportional_80", "families_detected_at_least_once_median"] == 250
    assert out.figure.exists()


def test_reciprocal_exclusion_definition():
    candidate = pd.DataFrame(
        [{"target": "p1", "ali_from": 10, "ali_to": 100, "domain_score": 50.0}]
    )
    canonical = pd.DataFrame(
        [{"target": "p1", "ali_from": 20, "ali_to": 90, "domain_score": 46.0}]
    )
    result = reciprocal_exclusion_grid(candidate, canonical, (0.5,), (5.0, 2.0))
    assert result.loc[result.score_tolerance_bits == 5.0, "retained"].item() == 0
    assert result.loc[result.score_tolerance_bits == 2.0, "retained"].item() == 1


def test_independent_projection_resampling(tmp_path):
    out = run_independent_projection_resampling(
        ROOT / "release/data/presence_absence_matrix.csv",
        ROOT / "release/tables/genome_accounting.csv",
        ROOT / "release/tables/discovery_panel.csv",
        tmp_path,
        replicates=5,
        seed=20_000,
    )
    summary = pd.read_csv(out.summary).iloc[0]
    assert summary["genomes_per_replicate"] == 331
    assert summary["discovery_projection_overlap"] == 0
    assert 230 <= summary["families_detected_median"] <= 250


def test_motif_helpers_preserve_length_and_rank():
    import numpy as np

    rng = np.random.default_rng(7)
    sequence = "ACGT" * 70
    shuffled = dinucleotide_shuffle(sequence, rng)
    assert len(shuffled) == len(sequence)
    assert auc(np.array([3.0, 4.0]), np.array([1.0, 2.0])) == 1.0


def test_quality_sensitivity_accounts_for_all_projection_genomes(tmp_path):
    out = run_quality_sensitivity(
        ROOT / "release/data/presence_absence_matrix.csv",
        ROOT / "release/data/projection_genome_metadata.csv",
        tmp_path,
    )
    summary = pd.read_csv(out.summary)
    origin = summary.loc[summary.stratification == "origin"]
    assert origin.n_genomes.sum() == 3192
    assert out.figure.exists()


def test_versioned_genus_blocked_results_have_no_genus_overlap():
    detail = pd.read_csv(ROOT / "results/genus_blocked_motif_replicates.csv")
    summary = pd.read_csv(ROOT / "results/genus_blocked_motif_summary.csv")
    assert len(detail) == 12 * 30
    assert detail.genus_overlap.eq(0).all()
    assert (summary.fraction_auc_above_0_5_vs_other_families >= (26 / 30)).sum() == 7


def test_taxonomic_clustered_prevalence_is_complete_and_bounded():
    summary = pd.read_csv(ROOT / "results/taxonomic_clustered_prevalence.csv")
    assert len(summary) == 250
    assert summary.n_genomes.eq(3192).all()
    assert summary.n_inferred_genera.eq(533).all()
    bounded = summary.filter(regex="prevalence|q025|q975")
    assert bounded.ge(0).all().all() and bounded.le(1).all().all()


def test_family_separability_is_explicitly_subset_limited():
    detail = pd.read_csv(ROOT / "results/family_separability_sequences.csv")
    summary = pd.read_csv(ROOT / "results/family_separability_summary.csv")
    assert len(detail) == 782
    assert len(summary) == 10
    assert summary.leave_one_out_accuracy.between(0, 1).all()


def test_lineage_blocked_motif_has_no_lineage_overlap():
    blocks = pd.read_csv(ROOT / "results/lineage_blocked_motif_blocks.csv")
    summary = pd.read_csv(ROOT / "results/lineage_blocked_motif_summary.csv")
    assert blocks.lineage_overlap.eq(0).all()
    assert (blocks.test_sequences >= 5).all()
    assert (blocks.train_sequences >= 20).all()
    assert (summary.evaluable_lineages > 0).sum() == 10


def test_extended_validation_registry_does_not_overclaim_completion():
    registry = pd.read_csv(ROOT / "release/tables/extended_validation_registry.csv")
    completed = set(registry.loc[registry.status == "complete", "item"])
    assert completed == {
        "Taxonomic clustered prevalence",
        "Ten-family sequence separability",
        "Genus-blocked motif validation",
        "Major-lineage-blocked motif validation",
    }
