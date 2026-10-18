from __future__ import annotations

import re
import sys
from pathlib import Path

from docx import Document


DOI = re.compile(r"https://doi\.org/(10\.\d{4,9}/[-._;()/:A-Za-z0-9]+)")


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("manuscript/Revised_archaeal_DBD_manuscript.docx")
    document = Document(path)
    references = [p.text.strip() for p in document.paragraphs if re.match(r"^\d+\. ", p.text.strip())]
    dois = []
    problems = []
    expected = 1
    for reference in references:
        number = int(reference.split(".", 1)[0])
        if number != expected:
            problems.append(f"Expected reference {expected}, found {number}")
        expected += 1
        match = DOI.search(reference)
        if not match:
            problems.append(f"Reference {number} has no DOI URL")
        else:
            dois.append(match.group(1).lower().rstrip("."))
    duplicates = sorted({doi for doi in dois if dois.count(doi) > 1})
    if duplicates:
        problems.append(f"Duplicate DOI values: {duplicates}")
    if len(references) != 40:
        problems.append(f"Expected 40 references, found {len(references)}")
    if problems:
        raise SystemExit("\n".join(problems))
    print(f"Reference audit passed: {len(references)} sequential entries, each with one unique DOI URL.")


if __name__ == "__main__":
    main()
