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
class LineageBlockedOutputs:
    blocks: Path
    summary: Path
    figure: Path


def run_lineage_blocked_motif(
    upstream_path: Path,
    motif_table_path: Path,
    discovery_panel_path: Path,
    output_dir: Path,
    seed: int = 91_007,
    minimum_test_sequences: int = 5,
) -> LineageBlockedOutputs:
    """Evaluate motif transfer by leaving one major lineage out at a time."""
    with gzip.open(upstream_path, "rt", encoding="utf-8") as handle:
        upstreams: dict[str, dict[str, str]] = json.load(handle)
    motifs = pd.read_csv(motif_table_path).set_index("family_id")
    panel = pd.read_csv(discovery_panel_path).set_index("accession")
    lineage_by_accession = panel["major_group"].astype(str).to_dict()
    records_by_family: dict[str, list[tuple[str, str, str]]] = {}
    for family, records in upstreams.items():
        records_by_family[family] = [
            (identifier, sequence.upper(), lineage_by_accession[identifier.split("|")[0]])
            for identifier, sequence in records.items()
            if len(sequence) >= 250 and identifier.split("|")[0] in lineage_by_accession
        ]

    rows = []
    for family, records in sorted(records_by_family.items()):
        width = len(str(motifs.at[family, "consensus"]))
        lineages = sorted({lineage for _, _, lineage in records})
        other_records = [record for other, values in records_by_family.items() if other != family for record in values]
        for block_index, heldout in enumerate(lineages):
            train = [sequence for _, sequence, lineage in records if lineage != heldout]
            test = [sequence for _, sequence, lineage in records if lineage == heldout]
            controls = [sequence for _, sequence, lineage in other_records if lineage == heldout]
            if len(test) < minimum_test_sequences or len(train) < 20 or len(controls) < len(test):
                continue
            rng = np.random.default_rng(seed + 1000 * int(family.split("-")[1]) + block_index)
            pwm, background = _fit_pwm(train, width, rng)
            positive = score_sequences(test, pwm, background)
            shuffled = score_sequences([dinucleotide_shuffle(sequence, rng) for sequence in test], pwm, background)
            target_gc = np.mean([(sequence.count("G") + sequence.count("C")) / len(sequence) for sequence in test])
            ranked = sorted(controls, key=lambda sequence: abs((sequence.count("G") + sequence.count("C")) / len(sequence) - target_gc))
            matched = ranked[: len(test)]
            other = score_sequences(matched, pwm, background)
            rows.append(
                {
                    "family_id": family,
                    "heldout_lineage": heldout,
                    "train_lineages": len(lineages) - 1,
                    "lineage_overlap": 0,
                    "train_sequences": len(train),
                    "test_sequences": len(test),
                    "control_sequences": len(matched),
                    "auc_vs_dinucleotide_shuffle": auc(positive, shuffled),
                    "auc_vs_gc_matched_other_families": auc(positive, other),
                }
            )
    detail = pd.DataFrame(rows)
    if detail.empty:
        raise ValueError("No lineage blocks met the preregistered sample-size requirements")
    summaries = []
    for family in sorted(records_by_family):
        subset = detail.loc[detail.family_id == family]
        summaries.append(
            {
                "family_id": family,
                "evaluable_lineages": len(subset),
                "total_test_sequences": int(subset["test_sequences"].sum()) if len(subset) else 0,
                "median_auc_vs_shuffle": float(subset["auc_vs_dinucleotide_shuffle"].median()) if len(subset) else np.nan,
                "median_auc_vs_other_families": float(subset["auc_vs_gc_matched_other_families"].median()) if len(subset) else np.nan,
                "lineages_auc_above_0_5_vs_shuffle": int((subset["auc_vs_dinucleotide_shuffle"] > 0.5).sum()),
                "lineages_auc_above_0_5_vs_other_families": int((subset["auc_vs_gc_matched_other_families"] > 0.5).sum()),
            }
        )
    summary = pd.DataFrame(summaries)
    output_dir.mkdir(parents=True, exist_ok=True)
    blocks_path = output_dir / "lineage_blocked_motif_blocks.csv"
    summary_path = output_dir / "lineage_blocked_motif_summary.csv"
    detail.to_csv(blocks_path, index=False)
    summary.to_csv(summary_path, index=False)

    figure_data = detail.sort_values(["family_id", "heldout_lineage"])
    fig, axis = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    for family, subset in figure_data.groupby("family_id"):
        axis.scatter(
            subset["auc_vs_dinucleotide_shuffle"],
            subset["auc_vs_gc_matched_other_families"],
            s=np.clip(subset["test_sequences"], 15, 120),
            alpha=0.7,
            label=family,
        )
    axis.axhline(0.5, linestyle="--", linewidth=1, color="#777777")
    axis.axvline(0.5, linestyle="--", linewidth=1, color="#777777")
    axis.set(xlim=(0, 1), ylim=(0, 1))
    axis.set_xlabel("AUC versus dinucleotide shuffle")
    axis.set_ylabel("AUC versus same-lineage other families")
    axis.set_title("Leave-one-major-lineage-out motif transfer")
    axis.legend(ncol=2, fontsize=7, frameon=False)
    figure_path = output_dir / "FigureS13_lineage_blocked_motif_validation.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return LineageBlockedOutputs(blocks_path, summary_path, figure_path)
