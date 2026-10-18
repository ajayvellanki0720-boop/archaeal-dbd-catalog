from __future__ import annotations

import gzip
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata


BASES = "ACGT"
BASE_INDEX = {base: index for index, base in enumerate(BASES)}


@dataclass(frozen=True)
class MotifHoldoutOutputs:
    replicates: Path
    summary: Path
    figure: Path


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _one_hot(sequences: list[str], length: int = 250) -> np.ndarray:
    kept = [sequence.upper()[-length:] for sequence in sequences if len(sequence) >= length]
    encoded = np.zeros((len(kept), length, 4), dtype=np.float32)
    for row, sequence in enumerate(kept):
        for column, base in enumerate(sequence):
            index = BASE_INDEX.get(base)
            if index is not None:
                encoded[row, column, index] = 1.0
    return encoded


def _fit_pwm(
    sequences: list[str], width: int, rng: np.random.Generator, starts: int = 3, iterations: int = 20
) -> tuple[np.ndarray, np.ndarray]:
    encoded = _one_hot(sequences)
    n_sequences, length, _ = encoded.shape
    positions = length - width + 1
    background_counts = encoded.sum(axis=(0, 1)) + 1.0
    background = background_counts / background_counts.sum()
    background_log = np.log(background)
    background_score = np.zeros((n_sequences, positions), dtype=np.float32)
    for offset in range(width):
        background_score += encoded[:, offset : offset + positions, :] @ background_log

    best_pwm = None
    best_likelihood = -np.inf
    for _ in range(starts):
        sequence_index = int(rng.integers(n_sequences))
        position = int(rng.integers(positions))
        pwm = np.full((width, 4), 0.1, dtype=np.float32)
        pwm += 0.6 * encoded[sequence_index, position : position + width, :]
        pwm /= pwm.sum(axis=1, keepdims=True)
        for _ in range(iterations):
            score = np.zeros((n_sequences, positions), dtype=np.float32)
            pwm_log = np.log(np.clip(pwm, 1e-8, None))
            for offset in range(width):
                score += encoded[:, offset : offset + positions, :] @ pwm_log[offset]
            likelihood_ratio = score - background_score
            likelihood_ratio -= likelihood_ratio.max(axis=1, keepdims=True)
            weights = np.exp(likelihood_ratio)
            weights /= weights.sum(axis=1, keepdims=True)
            counts = np.full((width, 4), 0.25, dtype=np.float32)
            for offset in range(width):
                counts[offset] += np.einsum(
                    "np,npb->b", weights, encoded[:, offset : offset + positions, :]
                )
            new_pwm = counts / counts.sum(axis=1, keepdims=True)
            if np.abs(new_pwm - pwm).max() < 1e-4:
                pwm = new_pwm
                break
            pwm = new_pwm
        final = _score_encoded(encoded, pwm, background)
        likelihood = float(final.sum())
        if likelihood > best_likelihood:
            best_likelihood = likelihood
            best_pwm = pwm.copy()
    assert best_pwm is not None
    return best_pwm, background


def _score_encoded(encoded: np.ndarray, pwm: np.ndarray, background: np.ndarray) -> np.ndarray:
    width = len(pwm)
    positions = encoded.shape[1] - width + 1
    pwm_log = np.log(np.clip(pwm, 1e-8, None))
    background_log = np.log(np.clip(background, 1e-8, None))
    forward = np.zeros((len(encoded), positions), dtype=np.float32)
    for offset in range(width):
        forward += encoded[:, offset : offset + positions, :] @ (pwm_log[offset] - background_log)
    reverse_pwm = pwm[::-1, ::-1]
    reverse_log = np.log(np.clip(reverse_pwm, 1e-8, None))
    reverse = np.zeros((len(encoded), positions), dtype=np.float32)
    for offset in range(width):
        reverse += encoded[:, offset : offset + positions, :] @ (reverse_log[offset] - background_log)
    return np.maximum(forward.max(axis=1), reverse.max(axis=1))


def score_sequences(sequences: list[str], pwm: np.ndarray, background: np.ndarray) -> np.ndarray:
    return _score_encoded(_one_hot(sequences), pwm, background)


def dinucleotide_shuffle(sequence: str, rng: np.random.Generator) -> str:
    sequence = "".join(base for base in sequence.upper() if base in BASE_INDEX)
    if len(sequence) < 4:
        return sequence
    edges: dict[str, list[str]] = defaultdict(list)
    for left, right in zip(sequence[:-1], sequence[1:]):
        edges[left].append(right)
    for _ in range(30):
        shuffled = {key: values[:] for key, values in edges.items()}
        for values in shuffled.values():
            rng.shuffle(values)
        result = [sequence[0]]
        current = sequence[0]
        for _ in range(len(sequence) - 1):
            if not shuffled[current]:
                break
            current = shuffled[current].pop()
            result.append(current)
        if len(result) == len(sequence):
            return "".join(result)
    values = list(sequence)
    rng.shuffle(values)
    return "".join(values)


def auc(positive: np.ndarray, negative: np.ndarray) -> float:
    values = np.concatenate([positive, negative])
    ranks = rankdata(values, method="average")
    n_positive = len(positive)
    n_negative = len(negative)
    rank_sum = ranks[:n_positive].sum()
    return float((rank_sum - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative))


def run_motif_holdout(
    upstream_path: Path,
    motif_table_path: Path,
    output_dir: Path,
    replicates: int = 30,
    seed: int = 30_000,
) -> MotifHoldoutOutputs:
    with gzip.open(upstream_path, "rt", encoding="utf-8") as handle:
        upstreams: dict[str, dict[str, str]] = json.load(handle)
    motif_table = pd.read_csv(motif_table_path).set_index("family_id")
    all_sequences = {
        family: [(identifier, sequence.upper()) for identifier, sequence in records.items() if len(sequence) >= 250]
        for family, records in upstreams.items()
    }

    rows = []
    for family, records in sorted(all_sequences.items()):
        width = len(str(motif_table.at[family, "consensus"]))
        genomes = sorted({identifier.split("|")[0] for identifier, _ in records})
        other_records = [
            record
            for other_family, other in all_sequences.items()
            if other_family != family
            for record in other
        ]
        for replicate in range(replicates):
            rng = np.random.default_rng(seed + replicate + 1_000 * int(family.split("-")[1]))
            shuffled_genomes = np.array(genomes, dtype=object)
            rng.shuffle(shuffled_genomes)
            test_count = max(1, int(round(0.3 * len(shuffled_genomes))))
            test_genomes = set(shuffled_genomes[:test_count])
            train = [sequence for identifier, sequence in records if identifier.split("|")[0] not in test_genomes]
            test = [sequence for identifier, sequence in records if identifier.split("|")[0] in test_genomes]
            pwm, background = _fit_pwm(train, width, rng)
            positive_scores = score_sequences(test, pwm, background)
            shuffled = [dinucleotide_shuffle(sequence, rng) for sequence in test]
            shuffled_scores = score_sequences(shuffled, pwm, background)

            target_gc = sum(sequence.count("G") + sequence.count("C") for sequence in test) / sum(map(len, test))
            gc_ranked = sorted(
                other_records,
                key=lambda item: abs(
                    (item[1].count("G") + item[1].count("C")) / len(item[1]) - target_gc
                ),
            )
            start = int(rng.integers(max(1, min(len(test), len(gc_ranked)))))
            matched = [sequence for _, sequence in gc_ranked[start : start + len(test)]]
            if len(matched) < len(test):
                matched.extend(sequence for _, sequence in gc_ranked[: len(test) - len(matched)])
            other_scores = score_sequences(matched, pwm, background)
            rows.append(
                {
                    "family_id": family,
                    "replicate": replicate,
                    "motif_width": width,
                    "train_genomes": len(genomes) - test_count,
                    "test_genomes": test_count,
                    "test_sequences": len(test),
                    "auc_vs_dinucleotide_shuffle": auc(positive_scores, shuffled_scores),
                    "auc_vs_gc_matched_other_families": auc(positive_scores, other_scores),
                    "median_score_difference_vs_shuffle": float(np.median(positive_scores) - np.median(shuffled_scores)),
                    "median_score_difference_vs_other": float(np.median(positive_scores) - np.median(other_scores)),
                }
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    detail = pd.DataFrame(rows)
    detail_path = output_dir / "motif_holdout_replicates.csv"
    detail.to_csv(detail_path, index=False)

    summary_rows = []
    for family, subset in detail.groupby("family_id", sort=True):
        successes_shuffle = int((subset.auc_vs_dinucleotide_shuffle > 0.5).sum())
        successes_other = int((subset.auc_vs_gc_matched_other_families > 0.5).sum())
        summary_rows.append(
            {
                "family_id": family,
                "replicates": len(subset),
                "median_auc_vs_dinucleotide_shuffle": float(subset.auc_vs_dinucleotide_shuffle.median()),
                "auc_vs_dinucleotide_shuffle_q025": float(subset.auc_vs_dinucleotide_shuffle.quantile(0.025)),
                "auc_vs_dinucleotide_shuffle_q975": float(subset.auc_vs_dinucleotide_shuffle.quantile(0.975)),
                "replicates_auc_above_0_5_vs_dinucleotide_shuffle": successes_shuffle,
                "fraction_auc_above_0_5_vs_dinucleotide_shuffle": successes_shuffle / len(subset),
                "median_auc_vs_gc_matched_other_families": float(subset.auc_vs_gc_matched_other_families.median()),
                "auc_vs_other_q025": float(subset.auc_vs_gc_matched_other_families.quantile(0.025)),
                "auc_vs_other_q975": float(subset.auc_vs_gc_matched_other_families.quantile(0.975)),
                "replicates_auc_above_0_5_vs_other_families": successes_other,
                "fraction_auc_above_0_5_vs_other_families": successes_other / len(subset),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary_path = output_dir / "motif_holdout_summary.csv"
    summary.to_csv(summary_path, index=False)

    order = summary.sort_values("median_auc_vs_gc_matched_other_families")["family_id"].tolist()
    positions = np.arange(len(order))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True, constrained_layout=True)
    for axis, metric, title in [
        (axes[0], "auc_vs_dinucleotide_shuffle", "Held-out sequence versus dinucleotide shuffle"),
        (axes[1], "auc_vs_gc_matched_other_families", "Held-out sequence versus other families"),
    ]:
        distributions = [detail.loc[detail.family_id == family, metric].to_numpy() for family in order]
        axis.boxplot(distributions, vert=False, tick_labels=order, showfliers=False)
        axis.axvline(0.5, color="#555555", linestyle="--", linewidth=1)
        axis.set_xlim(0, 1)
        axis.set_xlabel("Held-out AUC")
        axis.set_title(title)
    figure_path = output_dir / "FigureS7_motif_holdout_validation.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return MotifHoldoutOutputs(detail_path, summary_path, figure_path)
