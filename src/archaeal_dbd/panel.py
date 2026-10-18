from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCHEMES = ("proportional_80", "rare_full_major_80", "equal_cap_30")
LABELS = {
    "proportional_80": "80% within each lineage",
    "rare_full_major_80": "All rare lineages; 80% major caps",
    "equal_cap_30": "Maximum 30 per lineage",
}


@dataclass(frozen=True)
class PanelSensitivityOutputs:
    replicates: Path
    summary: Path
    figure: Path


def sample_panel(panel: pd.DataFrame, rng: np.random.Generator, scheme: str) -> list[str]:
    required = {"major_group", "accession"}
    missing = required - set(panel.columns)
    if missing:
        raise ValueError(f"Panel table is missing columns: {sorted(missing)}")
    chosen: list[str] = []
    for group, subset in panel.groupby("major_group", sort=True):
        if scheme == "proportional_80":
            n = max(1, int(round(0.8 * len(subset))))
        elif scheme == "equal_cap_30":
            n = min(30, len(subset))
        elif scheme == "rare_full_major_80":
            if group == "Euryarchaeota":
                n = min(112, len(subset))
            elif group == "Crenarchaeota (TACK)":
                n = min(56, len(subset))
            else:
                n = len(subset)
        else:
            raise ValueError(f"Unknown sampling scheme: {scheme}")
        indices = rng.choice(subset.index.to_numpy(), size=n, replace=False)
        chosen.extend(subset.loc[indices, "accession"].astype(str).tolist())
    return chosen


def run_panel_sensitivity(
    panel_path: Path,
    membership_path: Path,
    output_dir: Path,
    replicates: int = 500,
    seed: int = 10_000,
) -> PanelSensitivityOutputs:
    panel = pd.read_csv(panel_path)
    membership = pd.read_csv(membership_path)
    required = {"genome", "family_id"}
    missing = required - set(membership.columns)
    if missing:
        raise ValueError(f"Membership table is missing columns: {sorted(missing)}")
    families = membership["family_id"].dropna().unique()
    genome_family = membership[["genome", "family_id"]].drop_duplicates()

    rows: list[dict[str, float | int | str]] = []
    for scheme in SCHEMES:
        for replicate in range(replicates):
            rng = np.random.default_rng(seed + replicate)
            genomes = sample_panel(panel, rng, scheme)
            counts = (
                genome_family[genome_family["genome"].isin(genomes)]
                .groupby("family_id")["genome"]
                .nunique()
            )
            detected = int((counts >= 1).sum())
            retained = int((counts >= 5).sum())
            rows.append(
                {
                    "scheme": scheme,
                    "replicate": replicate,
                    "n_genomes": len(genomes),
                    "families_detected_at_least_once": detected,
                    "families_retaining_five_genomes": retained,
                    "fraction_detected_at_least_once": detected / len(families),
                    "fraction_retaining_five_genomes": retained / len(families),
                }
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    detail = pd.DataFrame(rows)
    detail_path = output_dir / "panel_resampling_sensitivity_replicates.csv"
    detail.to_csv(detail_path, index=False)

    summary_rows: list[dict[str, float | int | str]] = []
    for scheme, subset in detail.groupby("scheme", sort=True):
        row: dict[str, float | int | str] = {
            "scheme": scheme,
            "n_genomes": int(subset["n_genomes"].iloc[0]),
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
    summary_path = output_dir / "panel_resampling_sensitivity_summary.csv"
    summary.to_csv(summary_path, index=False)

    order = ["proportional_80", "rare_full_major_80", "equal_cap_30"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)
    specifications = [
        ("families_detected_at_least_once", "Families detected at least once"),
        ("families_retaining_five_genomes", "Families present in at least five genomes"),
    ]
    for axis, (metric, title) in zip(axes, specifications):
        values = [detail.loc[detail.scheme == scheme, metric].to_numpy() for scheme in order]
        axis.boxplot(values, tick_labels=[LABELS[s] for s in order], widths=0.55, showfliers=False)
        axis.axhline(len(families), color="#555555", linewidth=1, linestyle="--")
        axis.set_ylim(0, len(families) + 10)
        axis.set_ylabel(f"Families out of {len(families)}")
        axis.set_title(title)
        axis.tick_params(axis="x", rotation=20)
    fig.suptitle("Within-panel resampling sensitivity of the released family catalog")
    figure_path = output_dir / "FigureS5_panel_resampling_sensitivity.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return PanelSensitivityOutputs(detail_path, summary_path, figure_path)

