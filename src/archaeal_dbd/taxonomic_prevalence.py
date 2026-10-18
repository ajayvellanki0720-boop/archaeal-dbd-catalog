from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class TaxonomicPrevalenceOutputs:
    summary: Path
    figure: Path


def _genus_label(organism: object) -> str:
    text = str(organism).replace("[", "").replace("]", "").strip()
    tokens = text.split()
    if not tokens:
        return "unclassified"
    if tokens[0].lower() == "candidatus" and len(tokens) > 1:
        return f"Candidatus {tokens[1]}"
    return tokens[0]


def run_taxonomic_prevalence(
    presence_path: Path,
    metadata_path: Path,
    output_dir: Path,
    replicates: int = 1000,
    seed: int = 73_001,
) -> TaxonomicPrevalenceOutputs:
    """Estimate prevalence uncertainty by resampling inferred genera as clusters.

    Organism-name genera are a pragmatic surrogate because a versioned taxonomy
    table is not present. Results therefore quantify sensitivity to genome
    clustering, not fully phylogenetic comparative inference.
    """
    presence = pd.read_csv(presence_path, index_col=0)
    metadata = pd.read_csv(metadata_path).set_index("accession")
    if not set(presence.columns).issubset(metadata.index):
        raise ValueError("Projection metadata must cover every matrix accession")
    metadata = metadata.loc[presence.columns].copy()
    metadata["inferred_genus"] = metadata["organism"].map(_genus_label)
    genus_groups = {
        genus: np.flatnonzero(metadata["inferred_genus"].to_numpy() == genus)
        for genus in sorted(metadata["inferred_genus"].unique())
    }
    if len(genus_groups) < 10:
        raise ValueError("Too few inferred genera for clustered resampling")

    matrix = presence.to_numpy(dtype=float)
    genus_names = list(genus_groups)
    genus_means = np.column_stack([matrix[:, indices].mean(axis=1) for indices in genus_groups.values()])
    genus_sizes = np.array([len(indices) for indices in genus_groups.values()], dtype=float)
    rng = np.random.default_rng(seed)
    genus_balanced_boot = np.empty((replicates, matrix.shape[0]), dtype=float)
    cluster_weighted_boot = np.empty_like(genus_balanced_boot)
    for replicate in range(replicates):
        sampled = rng.integers(0, len(genus_names), size=len(genus_names))
        sampled_means = genus_means[:, sampled]
        sampled_sizes = genus_sizes[sampled]
        genus_balanced_boot[replicate] = sampled_means.mean(axis=1)
        cluster_weighted_boot[replicate] = np.average(sampled_means, axis=1, weights=sampled_sizes)

    rows = []
    for index, family in enumerate(presence.index):
        detected = genus_means[index] > 0
        rows.append(
            {
                "family_id": family,
                "n_genomes": matrix.shape[1],
                "n_inferred_genera": len(genus_names),
                "genera_with_family": int(detected.sum()),
                "genome_weighted_prevalence": float(matrix[index].mean()),
                "genus_balanced_prevalence": float(genus_means[index].mean()),
                "genus_balanced_q025": float(np.quantile(genus_balanced_boot[:, index], 0.025)),
                "genus_balanced_q975": float(np.quantile(genus_balanced_boot[:, index], 0.975)),
                "cluster_bootstrap_q025": float(np.quantile(cluster_weighted_boot[:, index], 0.025)),
                "cluster_bootstrap_q975": float(np.quantile(cluster_weighted_boot[:, index], 0.975)),
            }
        )
    summary = pd.DataFrame(rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "taxonomic_clustered_prevalence.csv"
    summary.to_csv(summary_path, index=False)

    fig, axis = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    axis.scatter(
        summary["genome_weighted_prevalence"],
        summary["genus_balanced_prevalence"],
        s=18,
        alpha=0.7,
        color="#2b6f9f",
    )
    axis.plot([0, 1], [0, 1], linestyle="--", linewidth=1, color="#555555")
    axis.set(xlim=(-0.02, 1.02), ylim=(-0.02, 1.02))
    axis.set_xlabel("Genome-weighted prevalence")
    axis.set_ylabel("Genus-balanced prevalence")
    axis.set_title("Sensitivity of family prevalence to genus imbalance")
    figure_path = output_dir / "FigureS11_taxonomic_clustered_prevalence.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return TaxonomicPrevalenceOutputs(summary_path, figure_path)
