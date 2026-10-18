from pathlib import Path

import pandas as pd

from archaeal_dbd.audit import run_release_audit
from archaeal_dbd.motif_holdout import auc, dinucleotide_shuffle
from archaeal_dbd.panel import run_panel_sensitivity
from archaeal_dbd.projection import run_independent_projection_resampling
from archaeal_dbd.thresholds import reciprocal_exclusion_grid


ROOT = Path(__file__).resolve().parents[1]


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
