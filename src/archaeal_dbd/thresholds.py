from __future__ import annotations

import pandas as pd


def reciprocal_exclusion_grid(
    candidates: pd.DataFrame,
    canonical: pd.DataFrame,
    overlap_thresholds=(0.3, 0.5, 0.7),
    score_tolerances=(0.0, 2.0, 5.0, 10.0, 20.0),
) -> pd.DataFrame:
    """Recompute candidate retention over overlap and bit-score thresholds.

    Expected columns are target, ali_from, ali_to, and domain_score. The function is
    included now, but the released package lacks the two raw hit tables needed to run it.
    """
    required = {"target", "ali_from", "ali_to", "domain_score"}
    for label, frame in (("candidates", candidates), ("canonical", canonical)):
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{label} table is missing columns: {sorted(missing)}")

    cand = candidates.sort_values("domain_score", ascending=False).drop_duplicates("target")
    canon = canonical.sort_values("domain_score", ascending=False).drop_duplicates("target")
    joined = cand.merge(canon, on="target", how="left", suffixes=("_candidate", "_canonical"))
    start = joined[["ali_from_candidate", "ali_from_canonical"]].max(axis=1)
    end = joined[["ali_to_candidate", "ali_to_canonical"]].min(axis=1)
    overlap = (end - start + 1).clip(lower=0)
    candidate_length = joined["ali_to_candidate"] - joined["ali_from_candidate"] + 1
    joined["overlap_fraction"] = overlap / candidate_length

    rows = []
    for threshold in overlap_thresholds:
        for tolerance in score_tolerances:
            explained = (
                joined["domain_score_canonical"].notna()
                & (joined["overlap_fraction"] > threshold)
                & (joined["domain_score_canonical"] >= joined["domain_score_candidate"] - tolerance)
            )
            rows.append(
                {
                    "overlap_threshold": threshold,
                    "score_tolerance_bits": tolerance,
                    "candidate_proteins": len(joined),
                    "excluded_by_canonical": int(explained.sum()),
                    "retained": int((~explained).sum()),
                    "retained_fraction": float((~explained).mean()),
                }
            )
    return pd.DataFrame(rows)

