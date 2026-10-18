from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .motif_holdout import _fit_pwm, auc, dinucleotide_shuffle, score_sequences


@dataclass(frozen=True)
class GenusBlockedOutputs:
    replicates: Path
    summary: Path
    figure: Path


def run_genus_blocked_motif(
    upstream_path: Path,
    motif_table_path: Path,
    discovery_panel_path: Path,
    output_dir: Path,
    replicates: int = 30,
    seed: int = 40_000,
) -> GenusBlockedOutputs:
    """Evaluate motif transfer to held-out genera without train/test genus overlap."""
    with gzip.open(upstream_path, "rt", encoding="utf-8") as handle:
        upstreams: dict[str, dict[str, str]] = json.load(handle)
    motif_table = pd.read_csv(motif_table_path).set_index("family_id")
    panel = pd.read_csv(discovery_panel_path).set_index("accession")
    accession_to_genus = panel["genus"].fillna("unclassified").astype(str).to_dict()

    records_by_family: dict[str, list[tuple[str, str, str]]] = {}
    for family, records in upstreams.items():
        cleaned = []
        for identifier, sequence in records.items():
            accession = identifier.split("|")[0]
            if len(sequence) >= 250 and accession in accession_to_genus:
                cleaned.append((identifier, sequence.upper(), accession_to_genus[accession]))
        records_by_family[family] = cleaned

    rows: list[dict[str, object]] = []
    for family, records in sorted(records_by_family.items()):
        width = len(str(motif_table.at[family, "consensus"]))
        genera = sorted({genus for _, _, genus in records})
        if len(genera) < 3:
            raise ValueError(f"{family} has fewer than three genera")
        other_records = [record for other, values in records_by_family.items() if other != family for record in values]
        family_number = int(family.split("-")[1])
        for replicate in range(replicates):
            rng = np.random.default_rng(seed + replicate + 1_000 * family_number)
            shuffled_genera = np.array(genera, dtype=object)
            rng.shuffle(shuffled_genera)
            test_count = max(1, int(round(0.3 * len(genera))))
            test_genera = set(shuffled_genera[:test_count])
            train = [sequence for _, sequence, genus in records if genus not in test_genera]
            test = [sequence for _, sequence, genus in records if genus in test_genera]
            pwm, background = _fit_pwm(train, width, rng)
            positive_scores = score_sequences(test, pwm, background)
            shuffled_scores = score_sequences([dinucleotide_shuffle(sequence, rng) for sequence in test], pwm, background)

            eligible_controls = [record for record in other_records if record[2] not in set(genera) - test_genera]
            if len(eligible_controls) < len(test):
                eligible_controls = other_records
            target_gc = sum(sequence.count("G") + sequence.count("C") for sequence in test) / sum(map(len, test))
            ranked = sorted(
                eligible_controls,
                key=lambda record: abs((record[1].count("G") + record[1].count("C")) / len(record[1]) - target_gc),
            )
            offset = int(rng.integers(max(1, min(len(test), len(ranked)))))
            matched = [sequence for _, sequence, _ in ranked[offset : offset + len(test)]]
            if len(matched) < len(test):
                matched.extend(sequence for _, sequence, _ in ranked[: len(test) - len(matched)])
            other_scores = score_sequences(matched, pwm, background)
            rows.append(
                {
                    "family_id": family,
                    "replicate": replicate,
                    "motif_width": width,
                    "train_genera": len(genera) - test_count,
                    "test_genera": test_count,
                    "genus_overlap": 0,
                    "test_sequences": len(test),
                    "auc_vs_dinucleotide_shuffle": auc(positive_scores, shuffled_scores),
                    "auc_vs_gc_matched_other_families": auc(positive_scores, other_scores),
                }
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    detail = pd.DataFrame(rows)
    detail_path = output_dir / "genus_blocked_motif_replicates.csv"
    detail.to_csv(detail_path, index=False)
    summary_rows = []
    for family, subset in detail.groupby("family_id", sort=True):
        row: dict[str, object] = {"family_id": family, "replicates": len(subset), "genus_overlap": 0}
        for metric, short in [
            ("auc_vs_dinucleotide_shuffle", "shuffle"),
            ("auc_vs_gc_matched_other_families", "other_families"),
        ]:
            successes = int((subset[metric] > 0.5).sum())
            row[f"median_auc_vs_{short}"] = float(subset[metric].median())
            row[f"auc_vs_{short}_q025"] = float(subset[metric].quantile(0.025))
            row[f"auc_vs_{short}_q975"] = float(subset[metric].quantile(0.975))
            row[f"replicates_auc_above_0_5_vs_{short}"] = successes
            row[f"fraction_auc_above_0_5_vs_{short}"] = successes / len(subset)
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary_path = output_dir / "genus_blocked_motif_summary.csv"
    summary.to_csv(summary_path, index=False)

    order = summary.sort_values("median_auc_vs_other_families")["family_id"].tolist()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True, constrained_layout=True)
    for axis, metric, title in [
        (axes[0], "auc_vs_dinucleotide_shuffle", "Held-out genera versus shuffled sequence"),
        (axes[1], "auc_vs_gc_matched_other_families", "Held-out genera versus other families"),
    ]:
        distributions = [detail.loc[detail.family_id == family, metric].to_numpy() for family in order]
        axis.boxplot(distributions, vert=False, tick_labels=order, showfliers=False)
        axis.axvline(0.5, color="#555555", linestyle="--", linewidth=1)
        axis.set_xlim(0, 1)
        axis.set_xlabel("Genus-blocked held-out AUC")
        axis.set_title(title)
    figure_path = output_dir / "FigureS9_genus_blocked_motif_validation.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return GenusBlockedOutputs(detail_path, summary_path, figure_path)
