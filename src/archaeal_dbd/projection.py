from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MATCHED_COUNTS = {
    "Euryarchaeota": 140,
    "Crenarchaeota_TACK": 70,
    "Thaumarchaeota_TACK": 51,
    "DPANN": 33,
    "Asgardarchaeota": 19,
    "Other_deep_TACK": 15,
    "Korarchaeota_TACK": 3,
}


@dataclass(frozen=True)
class ProjectionOutputs:
    replicates: Path
    summary: Path
    family_stability: Path
    figure: Path


def run_independent_projection_resampling(
    presence_path: Path,
    accounting_path: Path,
    discovery_panel_path: Path,
    output_dir: Path,
    replicates: int = 500,
    seed: int = 20_000,
) -> ProjectionOutputs:
    """Test catalog detection in lineage-matched projection genomes.

    This is independent-genome validation of the released family HMM projections,
    not de novo family discovery. Discovery accessions are explicitly excluded.
    """
    presence = pd.read_csv(presence_path, index_col=0)
    accounting = pd.read_csv(accounting_path).set_index("accession")
    discovery = set(pd.read_csv(discovery_panel_path)["accession"].astype(str))
    overlap = discovery.intersection(presence.columns)
    if overlap:
        presence = presence.drop(columns=sorted(overlap))

    available = [column for column in presence.columns if column in accounting.index]
    pools: dict[str, list[str]] = {}
    for group, requested in MATCHED_COUNTS.items():
        pool = [column for column in available if accounting.at[column, "group"] == group]
        if len(pool) < requested:
            raise ValueError(f"Only {len(pool)} {group} genomes are available; {requested} required")
        pools[group] = pool

    replicate_rows = []
    family_hits = np.zeros((len(presence), replicates), dtype=np.uint8)
    for replicate in range(replicates):
        rng = np.random.default_rng(seed + replicate)
        selected: list[str] = []
        for group, count in MATCHED_COUNTS.items():
            selected.extend(rng.choice(pools[group], count, replace=False).tolist())
        family_counts = presence[selected].sum(axis=1).to_numpy()
        family_hits[:, replicate] = (family_counts >= 1).astype(np.uint8)
        replicate_rows.append(
            {
                "replicate": replicate,
                "n_genomes": len(selected),
                "families_detected_at_least_once": int((family_counts >= 1).sum()),
                "families_present_in_five_genomes": int((family_counts >= 5).sum()),
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    detail = pd.DataFrame(replicate_rows)
    detail_path = output_dir / "independent_projection_resampling_replicates.csv"
    detail.to_csv(detail_path, index=False)

    summary = pd.DataFrame(
        [
            {
                "sampling_frame": "lineage-matched independent projection genomes",
                "replicates": replicates,
                "genomes_per_replicate": int(sum(MATCHED_COUNTS.values())),
                "discovery_projection_overlap": len(overlap),
                "families_detected_median": float(detail.families_detected_at_least_once.median()),
                "families_detected_q025": float(detail.families_detected_at_least_once.quantile(0.025)),
                "families_detected_q975": float(detail.families_detected_at_least_once.quantile(0.975)),
                "families_five_genomes_median": float(detail.families_present_in_five_genomes.median()),
                "families_five_genomes_q025": float(detail.families_present_in_five_genomes.quantile(0.025)),
                "families_five_genomes_q975": float(detail.families_present_in_five_genomes.quantile(0.975)),
            }
        ]
    )
    summary_path = output_dir / "independent_projection_resampling_summary.csv"
    summary.to_csv(summary_path, index=False)

    family = pd.DataFrame(
        {
            "family": presence.index,
            "detection_frequency": family_hits.mean(axis=1),
            "replicates_detected": family_hits.sum(axis=1),
        }
    ).sort_values(["detection_frequency", "family"])
    family_path = output_dir / "independent_projection_family_stability.csv"
    family.to_csv(family_path, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), constrained_layout=True)
    axes[0].hist(detail.families_detected_at_least_once, bins=np.arange(236.5, 251.5), color="#2878B5")
    axes[0].set_xlabel("Families detected per 331-genome replicate")
    axes[0].set_ylabel("Replicates")
    axes[0].set_title("Independent projection genomes")
    axes[1].hist(family.detection_frequency, bins=np.linspace(0, 1, 21), color="#41AB5D")
    axes[1].set_xlabel("Per-family detection frequency")
    axes[1].set_ylabel("Families")
    axes[1].set_title("Stability across 500 replicates")
    figure_path = output_dir / "FigureS6_independent_projection_resampling.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return ProjectionOutputs(detail_path, summary_path, family_path, figure_path)

