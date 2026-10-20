from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    path = ROOT / "release/tables/genome_accounting.csv"
    frame = pd.read_csv(path)
    hit_fates = {"analyzed_hit", "analyzed_hit_diamond"}
    frame["got_hit"] = frame["fate"].isin(hit_fates)
    frame.to_csv(path, index=False)


if __name__ == "__main__":
    main()
