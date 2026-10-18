from __future__ import annotations

from pathlib import Path

import pandas as pd

from .thresholds import reciprocal_exclusion_grid


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    example = root / "examples/minimal"
    candidate = pd.read_csv(example / "candidate_hits.csv")
    canonical = pd.read_csv(example / "canonical_hits.csv")
    result = reciprocal_exclusion_grid(candidate, canonical, (0.5,), (5.0, 2.0))
    expected = pd.read_csv(example / "expected_retention.csv")
    observed = result[["overlap_threshold", "score_tolerance_bits", "retained"]].reset_index(drop=True)
    if not observed.equals(expected):
        raise SystemExit("Smoke test failed: output differs from the checked-in expectation")
    print("Smoke test passed: reciprocal exclusion reproduced the expected two-row result.")


if __name__ == "__main__":
    main()
