from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import tempfile
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def protein_url(raw_url: str) -> str:
    """Repair the archived NCBI directory URL into a protein FASTA URL."""
    if "//_protein.faa.gz" not in raw_url:
        return raw_url
    directory = raw_url.split("//_protein.faa.gz", 1)[0].rstrip("/")
    assembly_name = directory.rsplit("/", 1)[-1]
    return f"{directory}/{assembly_name}_protein.faa.gz"


def download(url: str, destination: Path) -> tuple[str, int]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "archaeal-dbd-catalog/0.5.0 (sequence provenance recovery)"},
    )
    last_error: Exception | None = None
    for attempt in range(5):
        digest = hashlib.sha256()
        size = 0
        try:
            with urllib.request.urlopen(request, timeout=120) as response, destination.open("wb") as stream:
                while chunk := response.read(1024 * 1024):
                    stream.write(chunk)
                    digest.update(chunk)
                    size += len(chunk)
            return digest.hexdigest(), size
        except Exception as exc:
            last_error = exc
            destination.unlink(missing_ok=True)
            time.sleep(2**attempt)
    assert last_error is not None
    raise last_error


def fasta_records(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        identifier = None
        sequence: list[str] = []
        for line in stream:
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(sequence)
                identifier = line[1:].split(None, 1)[0]
                sequence = []
            else:
                sequence.append(line.strip())
        if identifier is not None:
            yield identifier, "".join(sequence)


def fasta_text_records(text: str):
    identifier = None
    sequence: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if identifier is not None:
                yield identifier, "".join(sequence)
            identifier = line[1:].split(None, 1)[0]
            sequence = []
        else:
            sequence.append(line.strip())
    if identifier is not None:
        yield identifier, "".join(sequence)


def fetch_proteins(accessions: list[str]) -> tuple[str, dict[str, str]]:
    query = urllib.parse.urlencode(
        {"db": "protein", "id": ",".join(accessions), "rettype": "fasta", "retmode": "text"}
    )
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "archaeal-dbd-catalog/0.5.0 (sequence provenance recovery)"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        text = response.read().decode("utf-8")
    return url, dict(fasta_text_records(text))


def write_deterministic_gzip(path: Path, text: str) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write(text.encode("utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recover released catalog member sequences from versioned NCBI assemblies"
    )
    parser.add_argument(
        "--panel", type=Path, default=ROOT / "release/tables/discovery_panel.csv"
    )
    parser.add_argument(
        "--membership", type=Path, default=ROOT / "release/data/novel_family_membership.csv"
    )
    parser.add_argument(
        "--output", type=Path, default=ROOT / "release/data/novel_family_members.faa.gz"
    )
    parser.add_argument(
        "--manifest", type=Path, default=ROOT / "release/tables/catalog_sequence_retrieval.csv"
    )
    parser.add_argument(
        "--missing", type=Path, default=ROOT / "release/tables/catalog_sequence_missing.csv"
    )
    parser.add_argument(
        "--fallback-manifest",
        type=Path,
        default=ROOT / "release/tables/catalog_sequence_fallback.csv",
    )
    parser.add_argument("--retrieval-date", default=date.today().isoformat())
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--cache-dir", type=Path)
    args = parser.parse_args()

    panel = pd.read_csv(args.panel)
    membership = pd.read_csv(args.membership)
    required_panel = {"accession", "ftp_protein_faa"}
    required_membership = {"family_id", "member", "genome"}
    if missing := required_panel - set(panel.columns):
        raise ValueError(f"Panel table is missing columns: {sorted(missing)}")
    if missing := required_membership - set(membership.columns):
        raise ValueError(f"Membership table is missing columns: {sorted(missing)}")

    membership = membership.copy()
    membership["protein"] = membership["member"].str.rsplit("|", n=1).str[-1]
    targets: dict[str, set[str]] = {
        genome: set(subset["protein"].astype(str))
        for genome, subset in membership.groupby("genome", sort=True)
    }
    urls = {
        str(row.accession): protein_url(str(row.ftp_protein_faa))
        for row in panel.itertuples(index=False)
        if str(row.accession) in targets
    }
    absent_genomes = sorted(set(targets) - set(urls))
    if absent_genomes:
        raise ValueError(f"Membership genomes absent from panel metadata: {absent_genomes[:5]}")

    temporary = None
    if args.cache_dir is None:
        temporary = tempfile.TemporaryDirectory(prefix="archaeal-dbd-proteomes-")
        cache_dir = Path(temporary.name)
    else:
        cache_dir = args.cache_dir
        cache_dir.mkdir(parents=True, exist_ok=True)

    downloads: dict[str, tuple[Path, str, int]] = {}
    failures: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        for genome, url in urls.items():
            path = cache_dir / f"{genome}_protein.faa.gz"
            if path.exists():
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                downloads[genome] = (path, digest, path.stat().st_size)
                continue
            futures[executor.submit(download, url, path)] = (genome, path)
        for future in as_completed(futures):
            genome, path = futures[future]
            try:
                digest, size = future.result()
                downloads[genome] = (path, digest, size)
            except Exception as exc:  # report every inaccessible versioned assembly
                failures[genome] = str(exc)

    recovered: dict[tuple[str, str], str] = {}
    manifest_rows = []
    missing_rows = []
    for genome in sorted(targets):
        wanted = targets[genome]
        found: dict[str, str] = {}
        if genome in downloads:
            path, digest, size = downloads[genome]
            for identifier, sequence in fasta_records(path):
                if identifier in wanted:
                    found[identifier] = sequence
        for protein in sorted(wanted):
            if protein in found:
                recovered[(genome, protein)] = found[protein]
            else:
                missing_rows.append(
                    {
                        "genome": genome,
                        "protein": protein,
                        "reason": failures.get(genome, "protein accession absent from assembly FASTA"),
                    }
                )
        manifest_rows.append(
            {
                "genome": genome,
                "source_url": urls[genome],
                "source_sha256": downloads.get(genome, (None, "", 0))[1],
                "source_bytes": downloads.get(genome, (None, "", 0))[2],
                "target_members": len(wanted),
                "recovered_members": len(found),
                "retrieval_date": args.retrieval_date,
                "status": "download_failed" if genome in failures else "complete" if len(found) == len(wanted) else "partial",
            }
        )

    fallback_rows = []
    initially_missing = {(row["genome"], row["protein"]) for row in missing_rows}
    missing_proteins = sorted({protein for _, protein in initially_missing})
    for start in range(0, len(missing_proteins), 100):
        batch = missing_proteins[start : start + 100]
        source_url, sequences = fetch_proteins(batch)
        for protein in batch:
            sequence = sequences.get(protein)
            fallback_rows.append(
                {
                    "protein": protein,
                    "source_url": source_url,
                    "sequence_sha256": hashlib.sha256(sequence.encode("utf-8")).hexdigest()
                    if sequence
                    else "",
                    "length_aa": len(sequence) if sequence else 0,
                    "retrieval_date": args.retrieval_date,
                    "status": "recovered" if sequence else "unavailable",
                }
            )
            if sequence:
                for genome, candidate in initially_missing:
                    if candidate == protein:
                        recovered[(genome, protein)] = sequence

    missing_rows = [
        row for row in missing_rows if (row["genome"], row["protein"]) not in recovered
    ]
    for row in manifest_rows:
        genome = str(row["genome"])
        fallback_count = sum(
            1 for candidate in targets[genome] if (genome, candidate) in initially_missing and (genome, candidate) in recovered
        )
        row["fallback_members"] = fallback_count
        row["total_recovered_members"] = sum(
            1 for candidate in targets[genome] if (genome, candidate) in recovered
        )
        if row["total_recovered_members"] == row["target_members"]:
            row["status"] = "complete_with_protein_fallback" if fallback_count else "complete"

    fasta = io.StringIO()
    for row in membership.sort_values(["family_id", "genome", "protein"]).itertuples(index=False):
        key = (str(row.genome), str(row.protein))
        if key not in recovered:
            continue
        fasta.write(f">{row.family_id}|{row.genome}|{row.protein}\n")
        sequence = recovered[key]
        for start in range(0, len(sequence), 80):
            fasta.write(sequence[start : start + 80] + "\n")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.missing.parent.mkdir(parents=True, exist_ok=True)
    args.fallback_manifest.parent.mkdir(parents=True, exist_ok=True)
    write_deterministic_gzip(args.output, fasta.getvalue())
    pd.DataFrame(manifest_rows).to_csv(args.manifest, index=False, quoting=csv.QUOTE_MINIMAL)
    pd.DataFrame(missing_rows, columns=["genome", "protein", "reason"]).to_csv(
        args.missing, index=False, quoting=csv.QUOTE_MINIMAL
    )
    pd.DataFrame(
        fallback_rows,
        columns=[
            "protein",
            "source_url",
            "sequence_sha256",
            "length_aa",
            "retrieval_date",
            "status",
        ],
    ).to_csv(args.fallback_manifest, index=False, quoting=csv.QUOTE_MINIMAL)

    print(f"Recovered {len(recovered):,} of {len(membership):,} catalog members")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.manifest}")
    print(f"Wrote {args.missing}")
    print(f"Wrote {args.fallback_manifest}")
    if temporary is not None:
        temporary.cleanup()
    if missing_rows:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
