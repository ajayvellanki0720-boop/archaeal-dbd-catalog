from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class QualitySensitivityOutputs:
    prevalence: Path
    summary: Path
    figure: Path


def _quality_class(metadata: pd.DataFrame) -> pd.Series:
    complete = pd.to_numeric(metadata["completeness"], errors="coerce")
    contamination = pd.to_numeric(metadata["contamination"], errors="coerce")
    reported = complete.notna() & contamination.notna()
    values = np.full(len(metadata), "not_reported", dtype=object)
    values[reported.to_numpy()] = "reported_below_threshold"
    values[(reported & (complete >= 90) & (contamination <= 5)).to_numpy()] = "reported_>=90%_complete_<=5%_contaminated"
    return pd.Series(values, index=metadata.index)


def run_quality_sensitivity(
    presence_path: Path,
    metadata_path: Path,
    output_dir: Path,
) -> QualitySensitivityOutputs:
    """Describe catalog prevalence across assembly/source-quality strata.

    These comparisons are sensitivity descriptions, not causal estimates: genome
    source, lineage, assembly level, and reported quality are confounded.
    """
    presence = pd.read_csv(presence_path, index_col=0)
    metadata = pd.read_csv(metadata_path).set_index("accession")
    common = [accession for accession in presence.columns if accession in metadata.index]
    if len(common) != presence.shape[1]:
        raise ValueError("Projection metadata must cover every presence/absence column")
    metadata = metadata.loc[common].copy()
    metadata["quality_class"] = _quality_class(metadata)
    metadata["assembly_class"] = metadata["assembly_level"].replace(
        {"Complete Genome": "complete_or_chromosome", "Chromosome": "complete_or_chromosome"}
    )

    rows: list[dict[str, object]] = []
    for stratification, column in [
        ("origin", "origin"),
        ("assembly_level", "assembly_class"),
        ("reported_quality", "quality_class"),
    ]:
        for stratum, accessions in metadata.groupby(column, dropna=False).groups.items():
            accessions = list(accessions)
            prevalence = presence[accessions].mean(axis=1)
            detections = presence[accessions].sum(axis=1)
            for family in presence.index:
                rows.append(
                    {
                        "stratification": stratification,
                        "stratum": str(stratum),
                        "family_id": family,
                        "n_genomes": len(accessions),
                        "genomes_with_family": int(detections.loc[family]),
                        "prevalence": float(prevalence.loc[family]),
                    }
                )
    detail = pd.DataFrame(rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    prevalence_path = output_dir / "quality_stratified_family_prevalence.csv"
    detail.to_csv(prevalence_path, index=False)

    summary = (
        detail.groupby(["stratification", "stratum"], as_index=False)
        .agg(
            n_genomes=("n_genomes", "first"),
            families_detected=("genomes_with_family", lambda values: int((values > 0).sum())),
            median_family_prevalence=("prevalence", "median"),
            q025_family_prevalence=("prevalence", lambda values: float(values.quantile(0.025))),
            q975_family_prevalence=("prevalence", lambda values: float(values.quantile(0.975))),
        )
    )
    summary_path = output_dir / "quality_stratified_summary.csv"
    summary.to_csv(summary_path, index=False)

    plots = [
        ("origin", ["isolate", "MAG", "SAG"], "Genome source"),
        ("assembly_level", ["complete_or_chromosome", "Scaffold", "Contig"], "Assembly level"),
        (
            "reported_quality",
            ["reported_>=90%_complete_<=5%_contaminated", "reported_below_threshold", "not_reported"],
            "Reported quality metadata",
        ),
    ]
    display_names = {
        "complete_or_chromosome": "complete/\nchromosome",
        "reported_>=90%_complete_<=5%_contaminated": "reported high quality",
        "reported_below_threshold": "reported below threshold",
        "not_reported": "not reported",
    }
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    for axis, (stratification, order, title) in zip(axes, plots):
        subset = detail.loc[detail.stratification == stratification]
        values = [subset.loc[subset.stratum == stratum, "prevalence"].to_numpy() for stratum in order]
        labels = [f"{display_names.get(stratum, stratum)}\n(n={int(subset.loc[subset.stratum == stratum, 'n_genomes'].iloc[0])})" for stratum in order]
        axis.boxplot(values, tick_labels=labels, showfliers=False)
        axis.set_title(title)
        axis.set_ylabel("Per-family prevalence")
        axis.tick_params(axis="x", labelrotation=25)
        axis.set_ylim(-0.02, 1.02)
    fig.suptitle("Sensitivity of catalogue prevalence to genome source and assembly quality")
    figure_path = output_dir / "FigureS8_assembly_quality_sensitivity.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return QualitySensitivityOutputs(prevalence_path, summary_path, figure_path)
