from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_FRACTIONS = tuple(np.round(np.arange(0.1, 1.01, 0.1), 1))


@dataclass(frozen=True)
class CatalogAccumulationOutputs:
    replicates: Path
    summary: Path
    figure: Path


def _stratified_sample(
    panel: pd.DataFrame,
    fraction: float,
    rng: np.random.Generator,
) -> list[str]:
    chosen: list[str] = []
    for _, subset in panel.groupby("major_group", sort=True):
        n = len(subset) if fraction >= 1 else max(1, int(round(fraction * len(subset))))
        indices = rng.choice(subset.index.to_numpy(), size=n, replace=False)
        chosen.extend(subset.loc[indices, "accession"].astype(str).tolist())
    return chosen


def run_catalog_accumulation(
    panel_path: Path,
    membership_path: Path,
    output_dir: Path,
    replicates: int = 500,
    seed: int = 70_000,
    fractions: tuple[float, ...] = DEFAULT_FRACTIONS,
) -> CatalogAccumulationOutputs:
    panel = pd.read_csv(panel_path)
    membership = pd.read_csv(membership_path)
    panel_required = {"major_group", "accession"}
    membership_required = {"genome", "family_id"}
    if missing := panel_required - set(panel.columns):
        raise ValueError(f"Panel table is missing columns: {sorted(missing)}")
    if missing := membership_required - set(membership.columns):
        raise ValueError(f"Membership table is missing columns: {sorted(missing)}")

    genome_family = membership[["genome", "family_id"]].drop_duplicates()
    n_families = int(membership["family_id"].nunique())
    rows: list[dict[str, float | int]] = []
    for fraction_index, fraction in enumerate(fractions):
        if not 0 < fraction <= 1:
            raise ValueError("Sampling fractions must be in the interval (0, 1].")
        n_runs = 1 if fraction >= 1 else replicates
        for replicate in range(n_runs):
            rng = np.random.default_rng(seed + fraction_index * replicates + replicate)
            genomes = _stratified_sample(panel, fraction, rng)
            counts = (
                genome_family[genome_family["genome"].isin(genomes)]
                .groupby("family_id")["genome"]
                .nunique()
            )
            rows.append(
                {
                    "sampling_fraction": fraction,
                    "replicate": replicate,
                    "n_genomes": len(genomes),
                    "families_detected_at_least_once": int((counts >= 1).sum()),
                    "families_retaining_five_genomes": int((counts >= 5).sum()),
                }
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    detail = pd.DataFrame(rows)
    detail_path = output_dir / "catalog_accumulation_replicates.csv"
    detail.to_csv(detail_path, index=False)

    summary_rows: list[dict[str, float | int]] = []
    for fraction, subset in detail.groupby("sampling_fraction", sort=True):
        row: dict[str, float | int] = {
            "sampling_fraction": fraction,
            "n_genomes_median": float(subset["n_genomes"].median()),
            "replicates": len(subset),
        }
        for metric in (
            "families_detected_at_least_once",
            "families_retaining_five_genomes",
        ):
            row[f"{metric}_median"] = float(subset[metric].median())
            row[f"{metric}_q025"] = float(subset[metric].quantile(0.025))
            row[f"{metric}_q975"] = float(subset[metric].quantile(0.975))
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary_path = output_dir / "catalog_accumulation_summary.csv"
    summary.to_csv(summary_path, index=False)

    fig, axis = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    colors = {
        "families_detected_at_least_once": "#1f77b4",
        "families_retaining_five_genomes": "#d55e00",
    }
    labels = {
        "families_detected_at_least_once": "Detected in at least one sampled genome",
        "families_retaining_five_genomes": "Present in at least five sampled genomes",
    }
    x = summary["n_genomes_median"].to_numpy()
    for metric in colors:
        median = summary[f"{metric}_median"].to_numpy()
        low = summary[f"{metric}_q025"].to_numpy()
        high = summary[f"{metric}_q975"].to_numpy()
        axis.plot(x, median, marker="o", linewidth=2, color=colors[metric], label=labels[metric])
        axis.fill_between(x, low, high, color=colors[metric], alpha=0.18, linewidth=0)
    axis.axhline(n_families, color="#555555", linewidth=1, linestyle="--")
    axis.set_xlabel("Genomes sampled with lineage stratification")
    axis.set_ylabel(f"Families out of {n_families}")
    axis.set_ylim(0, n_families + 10)
    axis.set_title("Accumulation of released families with panel size")
    axis.legend(frameon=False, loc="lower right")
    figure_path = output_dir / "FigureS14_catalog_accumulation.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return CatalogAccumulationOutputs(detail_path, summary_path, figure_path)
