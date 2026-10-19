from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def build_publication_figures(release_dir: Path, output_dir: Path) -> list[Path]:
    """Create claim-calibrated replacements for the three interpretive figures."""
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    tiers = pd.read_csv(release_dir / "tables/panel_regulator_tiers.csv")
    counts = tiers.groupby(["group", "tier"]).size().unstack(fill_value=0)
    counts = counts.div(counts.sum(axis=1), axis=0)
    tier_order = [column for column in ["known19", "characterized_beyond19", "uncharacterized"] if column in counts]
    labels = {
        "known19": "19 reference families",
        "characterized_beyond19": "other characterized domains",
        "uncharacterized_DUF": "DUF-annotated",
        "uncharacterized": "without characterized family assignment",
    }
    colors = ["#2166ac", "#67a9cf", "#d73027"]
    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    bottom = np.zeros(len(counts))
    for column, color in zip(tier_order, colors):
        ax.bar(counts.index, counts[column], bottom=bottom, label=labels[column], color=color)
        bottom += counts[column].to_numpy()
    ax.set_ylabel("Fraction of targeted DNA-binding-domain proteins")
    ax.set_title("Annotation composition of the targeted archaeal DBD panel")
    group_labels = {
        "Crenarchaeota (TACK)": "Crenarchaeota/TACK",
        "Thaumarchaeota (TACK)": "Thaumarchaeota/TACK",
        "Korarchaeota (TACK)": "Korarchaeota/TACK",
        "Other deep TACK": "Other deep TACK",
    }
    ax.set_xticks(range(len(counts.index)), [group_labels.get(value, value) for value in counts.index], rotation=42, ha="right", fontsize=8)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    path = output_dir / "Figure2_targeted_DBD_annotation_composition.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    prevalence = pd.read_csv(release_dir / "tables/corrected_prevalence.csv").set_index("family")
    shown = prevalence.mean(axis=1).sort_values(ascending=False).head(40).index
    matrix = prevalence.loc[shown]
    fig, ax = plt.subplots(figsize=(8.3, 9.6), constrained_layout=True)
    image = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(matrix.columns)), [name.replace("_", " ") for name in matrix.columns], rotation=45, ha="right")
    ax.set_yticks(range(len(matrix.index)), matrix.index, fontsize=7)
    ax.set_title("Cross-lineage prevalence of catalogued DBD families\nTop 40 families by mean corrected prevalence")
    ax.set_xlabel("Major archaeal group")
    ax.set_ylabel("Catalogued family")
    fig.colorbar(image, ax=ax, label="Fraction of group genomes with a confident detection")
    path = output_dir / "Figure4_cross_lineage_catalogue_prevalence.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    catalog = pd.read_csv(release_dir / "tables/novel_family_catalog.csv")
    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    categories = [
        ("beyond19_characterized", "Characterized beyond reference 19", "#2166ac"),
        ("uncharacterized", "Candidate noncanonical family", "#d73027"),
    ]
    for category, label, color in categories:
        subset = catalog.loc[catalog.novelty == category]
        ax.scatter(subset.n_genomes, subset.n_phyla, s=np.clip(subset.n_members, 20, 500) / 3, alpha=0.65, label=label, color=color, edgecolor="white", linewidth=0.3)
    ax.set_xscale("log")
    ax.set_xlabel("Genomes containing a catalogued family")
    ax.set_ylabel("Major archaeal groups represented")
    ax.set_title("Prevalence and taxonomic breadth of catalogued DBD families")
    ax.legend(frameon=False)
    path = output_dir / "Figure3_catalogued_family_prevalence_breadth.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    enrichment = pd.read_csv(release_dir / "tables/neighbor_enrichment_significant.csv")
    top = enrichment.sort_values("padj").head(80)
    family_order = top.groupby("family_id").padj.min().sort_values().head(25).index
    category_order = top.groupby("category").padj.min().sort_values().head(14).index
    matrix = top.pivot_table(index="family_id", columns="category", values="fold_enrich", aggfunc="max")
    matrix = matrix.reindex(index=family_order, columns=category_order)
    values = np.log2(matrix.to_numpy())
    fig, ax = plt.subplots(figsize=(9.5, 7.0), constrained_layout=True)
    image = ax.imshow(values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=np.nanpercentile(values, 95))
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=55, ha="right")
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    ax.set_title("Exploratory genomic-neighborhood associations")
    ax.set_xlabel("Neighbor annotation category")
    ax.set_ylabel("Catalogued DBD family")
    fig.colorbar(image, ax=ax, label="log2 fold enrichment")
    path = output_dir / "Figure4_genomic_neighborhood_associations.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    comparison = pd.read_csv(release_dir / "tables/uncharacterized_db_comparison.csv")
    counts = comparison["resolution"].value_counts()
    display = [
        ("CDD sequence match", int(counts.get("resolved_by_CDD", 0))),
        ("Structure-only match", int(counts.get("structure_only", 0))),
        ("Unresolved by resources tested", int(counts.get("genuinely_unresolved", 0))),
    ]
    fig, ax = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
    labels2, values2 = zip(*display)
    bars = ax.barh(labels2, values2, color=["#3b73b9", "#81a6d1", "#c84f56"])
    ax.invert_yaxis()
    ax.set_xlabel("Families among 79 lacking a characterized Pfam assignment")
    ax.set_title("Resolution using CDD and a curated structural reference set")
    for bar, value in zip(bars, values2):
        ax.text(value + 0.6, bar.get_y() + bar.get_height() / 2, str(value), va="center", weight="bold")
    path = output_dir / "Figure6_uncharacterized_family_resolution.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)
    return outputs


def build_workflow_figure(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 3.8), constrained_layout=True)
    ax.axis("off")
    steps = [
        ("Targeted archaeal\nDBD screen", "#d9edf7"),
        ("Reciprocal canonical\ndomain exclusion", "#dff0d8"),
        ("Sequence clustering\nand family HMMs", "#fcf8e3"),
        ("Independent-genome\nprojection", "#f2dede"),
        ("Neighborhood, motif,\nphylogeny and structure", "#e8daef"),
        ("Sensitivity and\nheld-out validation", "#d6eaf8"),
    ]
    x_positions = np.linspace(0.09, 0.91, len(steps))
    for index, ((label, color), x) in enumerate(zip(steps, x_positions)):
        ax.text(x, 0.52, label, ha="center", va="center", fontsize=10, transform=ax.transAxes,
                bbox={"boxstyle": "round,pad=0.55", "facecolor": color, "edgecolor": "#555555"})
        if index < len(steps) - 1:
            ax.annotate("", xy=(x_positions[index + 1] - 0.07, 0.52), xytext=(x + 0.07, 0.52),
                        xycoords=ax.transAxes, arrowprops={"arrowstyle": "->", "lw": 1.5, "color": "#555555"})
    ax.text(0.5, 0.9, "Analysis and validation workflow", ha="center", va="center", fontsize=15, weight="bold", transform=ax.transAxes)
    ax.text(0.5, 0.12, "Interpretation is limited to computationally supported candidate DNA-binding families; functional claims require experiment.",
            ha="center", va="center", fontsize=9, transform=ax.transAxes)
    path = output_dir / "FigureS10_analysis_workflow.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path
