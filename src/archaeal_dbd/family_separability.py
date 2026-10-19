from __future__ import annotations

import gzip
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import beta


@dataclass(frozen=True)
class FamilySeparabilityOutputs:
    sequences: Path
    summary: Path
    figure: Path


def _features(sequence: str, vocabulary: dict[str, int], k: int) -> np.ndarray:
    vector = np.zeros(len(vocabulary), dtype=float)
    counts = Counter(sequence[i : i + k] for i in range(max(0, len(sequence) - k + 1)))
    for token, count in counts.items():
        if token in vocabulary:
            vector[vocabulary[token]] = count
    norm = np.linalg.norm(vector)
    return vector / norm if norm else vector


def _binomial_interval(successes: int, total: int) -> tuple[float, float]:
    low = 0.0 if successes == 0 else float(beta.ppf(0.025, successes, total - successes + 1))
    high = 1.0 if successes == total else float(beta.ppf(0.975, successes + 1, total - successes))
    return low, high


def run_family_separability(
    sequence_path: Path,
    output_dir: Path,
    k: int = 2,
) -> FamilySeparabilityOutputs:
    """Test leave-one-sequence-out separability for released sequence subsets.

    This is an alignment-free composition check for ten released families, not
    validation of all 250 catalog boundaries.
    """
    with gzip.open(sequence_path, "rt", encoding="utf-8") as handle:
        family_sequences: dict[str, dict[str, str]] = json.load(handle)
    records = [
        (family, identifier, sequence.upper())
        for family, members in sorted(family_sequences.items())
        for identifier, sequence in members.items()
    ]
    alphabet = sorted({character for _, _, sequence in records for character in sequence})
    vocabulary = {
        a + b: index
        for index, (a, b) in enumerate((a, b) for a in alphabet for b in alphabet)
    }
    if k != 2:
        raise ValueError("The versioned analysis uses dipeptide composition (k=2)")
    matrix = np.vstack([_features(sequence, vocabulary, k) for _, _, sequence in records])
    labels = np.array([family for family, _, _ in records], dtype=object)
    families = sorted(set(labels))
    family_indices = {family: np.flatnonzero(labels == family) for family in families}
    family_sums = {family: matrix[indices].sum(axis=0) for family, indices in family_indices.items()}

    rows = []
    for row_index, (family, identifier, sequence) in enumerate(records):
        scores: dict[str, float] = {}
        for candidate in families:
            centroid_sum = family_sums[candidate].copy()
            denominator = len(family_indices[candidate])
            if candidate == family:
                centroid_sum -= matrix[row_index]
                denominator -= 1
            centroid = centroid_sum / denominator
            norm = np.linalg.norm(centroid)
            scores[candidate] = float(matrix[row_index] @ (centroid / norm)) if norm else 0.0
        best = max(scores, key=scores.get)
        best_other = max(score for candidate, score in scores.items() if candidate != family)
        rows.append(
            {
                "family_id": family,
                "sequence_id": identifier,
                "length_aa": len(sequence),
                "assigned_family": best,
                "correct": best == family,
                "own_centroid_similarity": scores[family],
                "best_other_centroid_similarity": best_other,
                "separation_margin": scores[family] - best_other,
            }
        )
    detail = pd.DataFrame(rows)
    summaries = []
    for family, subset in detail.groupby("family_id", sort=True):
        correct = int(subset["correct"].sum())
        low, high = _binomial_interval(correct, len(subset))
        summaries.append(
            {
                "family_id": family,
                "n_sequences": len(subset),
                "leave_one_out_accuracy": correct / len(subset),
                "accuracy_ci95_low": low,
                "accuracy_ci95_high": high,
                "median_own_similarity": float(subset["own_centroid_similarity"].median()),
                "median_best_other_similarity": float(subset["best_other_centroid_similarity"].median()),
                "median_separation_margin": float(subset["separation_margin"].median()),
                "fraction_positive_margin": float((subset["separation_margin"] > 0).mean()),
            }
        )
    summary = pd.DataFrame(summaries)
    output_dir.mkdir(parents=True, exist_ok=True)
    detail_path = output_dir / "family_separability_sequences.csv"
    summary_path = output_dir / "family_separability_summary.csv"
    detail.to_csv(detail_path, index=False)
    summary.to_csv(summary_path, index=False)

    order = summary.sort_values("median_separation_margin")["family_id"].tolist()
    distributions = [detail.loc[detail.family_id == family, "separation_margin"] for family in order]
    fig, axis = plt.subplots(figsize=(7, 5), constrained_layout=True)
    axis.boxplot(distributions, vert=False, tick_labels=order, showfliers=False)
    axis.axvline(0, linestyle="--", linewidth=1, color="#555555")
    axis.set_xlabel("Own-family minus nearest other-family centroid similarity")
    axis.set_title("Sequence-family separability in the prespecified 10-family subset")
    figure_path = output_dir / "FigureS12_family_sequence_separability.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return FamilySeparabilityOutputs(detail_path, summary_path, figure_path)
