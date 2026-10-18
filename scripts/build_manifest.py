from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_PARTS = {".git", ".pytest_cache", "__pycache__", ".venv", "work"}
EXCLUDED_FILES = {"MANIFEST.sha256"}
EXCLUDED_PREFIXES = {("results", "generated")}


def main() -> None:
    rows = []
    for path in sorted(ROOT.rglob("*")):
        relative = path.relative_to(ROOT)
        if (
            not path.is_file()
            or path.name in EXCLUDED_FILES
            or EXCLUDED_PARTS.intersection(relative.parts)
            or any(relative.parts[: len(prefix)] == prefix for prefix in EXCLUDED_PREFIXES)
        ):
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(f"{digest}  {relative.as_posix()}")
    (ROOT / "MANIFEST.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
