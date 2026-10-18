from __future__ import annotations

import argparse
import json
from pathlib import Path


REQUIRED = {
    "raw_discovery": [
        "external/raw_discovery/candidate_pool.faa.gz",
        "external/raw_discovery/candidate_hmm_hits.tsv.gz",
        "external/raw_discovery/canonical_hmm_hits.tsv.gz",
        "external/raw_discovery/genome_manifest.csv",
    ],
    "corrected_projection": [
        "external/corrected_projection/presence_absence_matrix_3526.csv.gz",
        "external/corrected_projection/copynumber_matrix_3526.csv.gz",
    ],
    "independent_discovery": [
        "external/independent_discovery/genome_manifest.csv",
        "external/independent_discovery/proteome_checksums.sha256",
    ],
    "external_annotation": [
        "external/annotation/query_sequences.faa",
        "external/annotation/database_manifest.csv",
    ],
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Inventory inputs required for blocked validation stages")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--strict", action="store_true", help="Return nonzero when any input is absent")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"stages": {}, "ready": True}
    for stage, relative_paths in REQUIRED.items():
        missing = [path for path in relative_paths if not (args.root / path).is_file()]
        report["stages"][stage] = {"ready": not missing, "missing": missing}
        report["ready"] = report["ready"] and not missing
    rendered = json.dumps(report, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    if args.strict and not report["ready"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
