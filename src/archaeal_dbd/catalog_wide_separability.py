from __future__ import annotations

import gzip
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
class CatalogWideSeparabilityOutputs:
    sequences: Path
    families: Path
    summary: Path
    figure: Path


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        identifier: str | None = None
        sequence: list[str] = []
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(sequence).upper()
                full_identifier = line[1:].split()[0]
                identifier = full_identifier.partition("|")[2] or full_identifier
                sequence = []
            elif line:
                sequence.append(line)
        if identifier is not None:
            records[identifier] = "".join(sequence).upper()
    return records


def _features(sequence: str, vocabulary: dict[str, int]) -> np.ndarray:
    vector = np.zeros(len(vocabulary), dtype=np.float32)
    for token, count in Counter(sequence[i : i + 2] for i in range(len(sequence) - 1)).items():
        index = vocabulary.get(token)
        if index is not None:
            vector[index] = count
    norm = np.linalg.norm(vector)
    return vector / norm if norm else vector


def _interval(successes: int, total: int) -> tuple[float, float]:
    low = 0.0 if successes == 0 else float(beta.ppf(0.025, successes, total - successes + 1))
    high = 1.0 if successes == total else float(beta.ppf(0.975, successes + 1, total - successes))
    return low, high


def run_catalog_wide_separability(
    sequence_path: Path,
    membership_path: Path,
    panel_path: Path,
    output_dir: Path,
    min_training_sequences: int = 2,
) -> CatalogWideSeparabilityOutputs:
    """Evaluate all released families with complete-genus holdout.

    Each test sequence is compared with normalized dipeptide-composition
    centroids fitted after removing every sequence from the test sequence's
    inferred genus. This is a taxonomically blocked boundary check, not an HMM
    calibration or a substitute for de novo family discovery.
    """
    membership = pd.read_csv(membership_path)
    panel = pd.read_csv(panel_path)
    genus_by_genome = panel.set_index("accession")["genus"].fillna("unresolved").to_dict()
    sequences = _read_fasta(sequence_path)
    membership["genus"] = membership["genome"].map(genus_by_genome).fillna("unresolved")
    membership["sequence"] = membership["member"].map(sequences)
    if membership["sequence"].isna().any():
        missing = int(membership["sequence"].isna().sum())
        raise ValueError(f"{missing} catalog members are missing from the sequence archive")

    alphabet = sorted({residue for sequence in membership["sequence"] for residue in sequence})
    vocabulary = {a + b: i for i, (a, b) in enumerate((a, b) for a in alphabet for b in alphabet)}
    matrix = np.vstack([_features(sequence, vocabulary) for sequence in membership["sequence"]])
    family_names = sorted(membership["family_id"].unique(), key=lambda value: int(value.split("-")[1]))
    family_index = {family: index for index, family in enumerate(family_names)}
    family_codes = membership["family_id"].map(family_index).to_numpy()

    total_sums = np.zeros((len(family_names), matrix.shape[1]), dtype=np.float64)
    total_counts = np.zeros(len(family_names), dtype=int)
    np.add.at(total_sums, family_codes, matrix)
    np.add.at(total_counts, family_codes, 1)

    detail_rows: list[dict[str, object]] = []
    for genus, test_frame in membership.groupby("genus", sort=True):
        test_indices = test_frame.index.to_numpy()
        heldout_sums = np.zeros_like(total_sums)
        heldout_counts = np.zeros_like(total_counts)
        np.add.at(heldout_sums, family_codes[test_indices], matrix[test_indices])
        np.add.at(heldout_counts, family_codes[test_indices], 1)
        train_sums = total_sums - heldout_sums
        train_counts = total_counts - heldout_counts
        eligible = train_counts >= min_training_sequences
        centroids = np.zeros_like(train_sums)
        centroids[eligible] = train_sums[eligible] / train_counts[eligible, None]
        norms = np.linalg.norm(centroids, axis=1)
        nonzero = eligible & (norms > 0)
        centroids[nonzero] /= norms[nonzero, None]
        candidate_indices = np.flatnonzero(nonzero)
        if len(candidate_indices) < 2:
            continue
        test_matrix = matrix[test_indices].astype(np.float64, copy=False)
        candidate_centroids = centroids[candidate_indices]
        if not np.isfinite(test_matrix).all() or not np.isfinite(candidate_centroids).all():
            raise ValueError("Non-finite values encountered in catalog-wide feature matrices")
        # Some Accelerate/BLAS builds emit spurious floating-point warnings for
        # this finite cosine-similarity product; the explicit finite checks and
        # postcondition guard the calculation itself.
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            scores = test_matrix @ candidate_centroids.T
        if not np.isfinite(scores).all():
            raise ValueError("Non-finite cosine similarities encountered")
        for local_row, global_row in enumerate(test_indices):
            true_code = family_codes[global_row]
            if not nonzero[true_code]:
                continue
            best_position = int(np.argmax(scores[local_row]))
            assigned_code = int(candidate_indices[best_position])
            true_position = int(np.flatnonzero(candidate_indices == true_code)[0])
            true_score = float(scores[local_row, true_position])
            other_scores = scores[local_row].copy()
            other_scores[true_position] = -np.inf
            best_other = float(np.max(other_scores))
            detail_rows.append(
                {
                    "family_id": family_names[true_code],
                    "sequence_id": membership.at[global_row, "member"],
                    "genome": membership.at[global_row, "genome"],
                    "held_out_genus": genus,
                    "assigned_family": family_names[assigned_code],
                    "correct": assigned_code == true_code,
                    "own_centroid_similarity": true_score,
                    "best_other_centroid_similarity": best_other,
                    "separation_margin": true_score - best_other,
                    "own_training_sequences": int(train_counts[true_code]),
                    "candidate_families": int(len(candidate_indices)),
                }
            )

    detail = pd.DataFrame(detail_rows)
    family_rows: list[dict[str, object]] = []
    for family in family_names:
        subset = detail.loc[detail["family_id"] == family]
        total_members = int((membership["family_id"] == family).sum())
        if subset.empty:
            family_rows.append(
                {
                    "family_id": family,
                    "catalog_sequences": total_members,
                    "evaluable_sequences": 0,
                    "evaluable_genera": 0,
                    "accuracy": np.nan,
                    "accuracy_ci95_low": np.nan,
                    "accuracy_ci95_high": np.nan,
                    "median_separation_margin": np.nan,
                    "fraction_positive_margin": np.nan,
                }
            )
            continue
        correct = int(subset["correct"].sum())
        low, high = _interval(correct, len(subset))
        family_rows.append(
            {
                "family_id": family,
                "catalog_sequences": total_members,
                "evaluable_sequences": len(subset),
                "evaluable_genera": subset["held_out_genus"].nunique(),
                "accuracy": correct / len(subset),
                "accuracy_ci95_low": low,
                "accuracy_ci95_high": high,
                "median_separation_margin": float(subset["separation_margin"].median()),
                "fraction_positive_margin": float((subset["separation_margin"] > 0).mean()),
            }
        )
    families = pd.DataFrame(family_rows)
    evaluable = families.dropna(subset=["accuracy"])
    summary = pd.DataFrame(
        [
            {
                "catalog_families": len(family_names),
                "catalog_sequences": len(membership),
                "evaluable_families": len(evaluable),
                "evaluable_sequences": len(detail),
                "evaluable_genera": detail["held_out_genus"].nunique(),
                "micro_accuracy": float(detail["correct"].mean()),
                "macro_family_accuracy": float(evaluable["accuracy"].mean()),
                "median_family_accuracy": float(evaluable["accuracy"].median()),
                "families_accuracy_at_least_0_8": int((evaluable["accuracy"] >= 0.8).sum()),
                "families_positive_median_margin": int((evaluable["median_separation_margin"] > 0).sum()),
                "method": "leave-one-inferred-genus-out normalized dipeptide centroid",
                "interpretation": "taxonomically blocked composition check; not HMM calibration",
            }
        ]
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    detail_path = output_dir / "catalog_wide_genus_blocked_sequences.csv"
    family_path = output_dir / "catalog_wide_genus_blocked_families.csv"
    summary_path = output_dir / "catalog_wide_genus_blocked_summary.csv"
    figure_path = output_dir / "FigureS15_catalog_wide_genus_blocked_separability.png"
    detail.to_csv(detail_path, index=False)
    families.to_csv(family_path, index=False)
    summary.to_csv(summary_path, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), constrained_layout=True)
    axes[0].hist(evaluable["accuracy"], bins=np.linspace(0, 1, 11), color="#4472C4", edgecolor="white")
    axes[0].axvline(evaluable["accuracy"].median(), color="#333333", linestyle="--", linewidth=1)
    axes[0].set_xlabel("Family accuracy")
    axes[0].set_ylabel("Number of families")
    axes[0].set_title("Complete-genus holdout")
    axes[1].scatter(
        families["catalog_sequences"], families["accuracy"],
        s=np.clip(families["evaluable_genera"].fillna(0), 2, 40), alpha=0.65, color="#ED7D31",
    )
    axes[1].set_xscale("log")
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].set_xlabel("Catalog sequences per family")
    axes[1].set_ylabel("Family accuracy")
    axes[1].set_title("Accuracy by family size")
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return CatalogWideSeparabilityOutputs(detail_path, family_path, summary_path, figure_path)
