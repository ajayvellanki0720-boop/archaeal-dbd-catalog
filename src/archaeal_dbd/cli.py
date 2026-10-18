from __future__ import annotations

import argparse
import json
from pathlib import Path

from .audit import run_release_audit
from .motif_holdout import run_motif_holdout
from .panel import run_panel_sensitivity
from .projection import run_independent_projection_resampling


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify and analyze the archaeal DBD release")
    parser.add_argument(
        "command",
        choices=("audit", "panel-sensitivity", "independent-projection", "motif-holdout", "all"),
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


if __name__ == "__main__":
    main()
