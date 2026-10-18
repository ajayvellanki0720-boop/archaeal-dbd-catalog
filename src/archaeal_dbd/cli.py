from __future__ import annotations

import argparse
import json
from pathlib import Path

from .audit import run_release_audit
from .panel import run_panel_sensitivity


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify and analyze the archaeal DBD release")
    parser.add_argument("command", choices=("audit", "panel-sensitivity", "all"))
    parser.add_argument("--release-dir", type=Path, default=Path("release"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/generated"))
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=10_000)
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


if __name__ == "__main__":
    main()

