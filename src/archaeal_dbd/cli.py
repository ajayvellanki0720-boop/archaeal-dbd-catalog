from __future__ import annotations

import argparse
import json
from pathlib import Path

from .audit import run_release_audit
from .catalog_accumulation import run_catalog_accumulation
from .catalog_wide_separability import run_catalog_wide_separability
from .family_separability import run_family_separability
from .genus_blocked_motif import run_genus_blocked_motif
from .lineage_blocked_motif import run_lineage_blocked_motif
from .motif_holdout import run_motif_holdout
from .panel import run_panel_sensitivity
from .projection import run_independent_projection_resampling
from .publication_figures import build_publication_figures, build_workflow_figure
from .quality_sensitivity import run_quality_sensitivity
from .taxonomic_prevalence import run_taxonomic_prevalence


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify and analyze the archaeal DBD release")
    parser.add_argument(
        "command",
        choices=(
            "audit", "panel-sensitivity", "catalog-accumulation", "independent-projection", "motif-holdout",
            "genus-blocked-motif", "lineage-blocked-motif", "quality-sensitivity",
            "taxonomic-prevalence", "family-separability", "catalog-wide-separability",
            "publication-figures", "all",
        ),
    )
    parser.add_argument("--release-dir", type=Path, default=Path("release"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/generated"))
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=10_000)
    parser.add_argument("--motif-replicates", type=int, default=30)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.command in {"audit", "all"}:
        result = run_release_audit(args.release_dir, args.output_dir / "release_audit.json")
        print(json.dumps(result, indent=2))
    if args.command in {"panel-sensitivity", "all"}:
        outputs = run_panel_sensitivity(
            args.release_dir / "tables/discovery_panel.csv",
            args.release_dir / "data/novel_family_membership.csv",
            args.output_dir,
            replicates=args.replicates,
            seed=args.seed,
        )
        print(f"Wrote {outputs.summary}")
        print(f"Wrote {outputs.figure}")
    if args.command in {"catalog-accumulation", "all"}:
        outputs = run_catalog_accumulation(
            args.release_dir / "tables/discovery_panel.csv",
            args.release_dir / "data/novel_family_membership.csv",
            args.output_dir,
            replicates=args.replicates,
            seed=args.seed + 60_000,
        )
        print(f"Wrote {outputs.summary}")
        print(f"Wrote {outputs.figure}")
    if args.command in {"independent-projection", "all"}:
        outputs = run_independent_projection_resampling(
            args.release_dir / "data/presence_absence_matrix.csv",
            args.release_dir / "tables/genome_accounting.csv",
            args.release_dir / "tables/discovery_panel.csv",
            args.output_dir,
            replicates=args.replicates,
            seed=args.seed + 10_000,
        )
        print(f"Wrote {outputs.summary}")
        print(f"Wrote {outputs.figure}")
    if args.command in {"motif-holdout", "all"}:
        outputs = run_motif_holdout(
            args.release_dir / "data/family_upstreams_250bp.json.gz",
            args.release_dir / "tables/motif_null_analysis.csv",
            args.output_dir,
            replicates=args.motif_replicates,
            seed=args.seed + 20_000,
        )
        print(f"Wrote {outputs.summary}")
        print(f"Wrote {outputs.figure}")
    if args.command in {"genus-blocked-motif", "all"}:
        outputs = run_genus_blocked_motif(
            args.release_dir / "data/family_upstreams_250bp.json.gz",
            args.release_dir / "tables/motif_null_analysis.csv",
            args.release_dir / "tables/discovery_panel.csv",
            args.output_dir,
            replicates=args.motif_replicates,
            seed=args.seed + 30_000,
        )
        print(f"Wrote {outputs.summary}")
        print(f"Wrote {outputs.figure}")
    if args.command in {"quality-sensitivity", "all"}:
        outputs = run_quality_sensitivity(
            args.release_dir / "data/presence_absence_matrix.csv",
            args.release_dir / "data/projection_genome_metadata.csv",
            args.output_dir,
        )
        print(f"Wrote {outputs.summary}")
        print(f"Wrote {outputs.figure}")
    if args.command in {"taxonomic-prevalence", "all"}:
        outputs = run_taxonomic_prevalence(
            args.release_dir / "data/presence_absence_matrix.csv",
            args.release_dir / "data/projection_genome_metadata.csv",
            args.output_dir,
        )
        print(f"Wrote {outputs.summary}")
        print(f"Wrote {outputs.figure}")
    if args.command in {"family-separability", "all"}:
        outputs = run_family_separability(
            args.release_dir / "data/tree_family_seqs.json.gz",
            args.output_dir,
        )
        print(f"Wrote {outputs.summary}")
        print(f"Wrote {outputs.figure}")
    if args.command in {"catalog-wide-separability", "all"}:
        outputs = run_catalog_wide_separability(
            args.release_dir / "data/novel_family_members.faa.gz",
            args.release_dir / "data/novel_family_membership.csv",
            args.release_dir / "tables/discovery_panel.csv",
            args.output_dir,
        )
        print(f"Wrote {outputs.summary}")
        print(f"Wrote {outputs.figure}")
    if args.command in {"lineage-blocked-motif", "all"}:
        outputs = run_lineage_blocked_motif(
            args.release_dir / "data/family_upstreams_250bp.json.gz",
            args.release_dir / "tables/motif_null_analysis.csv",
            args.release_dir / "tables/discovery_panel.csv",
            args.output_dir,
        )
        print(f"Wrote {outputs.summary}")
        print(f"Wrote {outputs.figure}")
    if args.command in {"publication-figures", "all"}:
        for path in build_publication_figures(args.release_dir, args.output_dir):
            print(f"Wrote {path}")
        print(f"Wrote {build_workflow_figure(args.output_dir)}")


if __name__ == "__main__":
    main()
