# AI_CODEX: Generated the parallel FASTA extraction stage and boundary metadata.
"""Split one fixed-width FASTA record into independently readable chunks.

Each worker writes a sequence containing its non-overlapping core plus up to
200 context bases on both sides. Coordinates in metadata are 1-based inclusive,
while offsets into the chunk sequence are 0-based half-open.

Examples:
    python -m module_2_assignment.fasta_chunks index INPUT OUTPUT --chunks 50
    python -m module_2_assignment.fasta_chunks worker OUTPUT 0
    python -m module_2_assignment.fasta_chunks local INPUT OUTPUT \
        --chunks 50 --workers 10
"""

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

CONTEXT_BASES = 200


def _write_json(path: Path, value: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")
    temporary.replace(path)


def index_fasta(fasta_path: str | Path, output_dir: str | Path, chunks: int) -> dict:
    """Validate a single-record, fixed-width FASTA and save seek information."""
    if not 1 <= chunks <= 50:
        raise ValueError("chunks must be between 1 and 50")
    fasta_path = Path(fasta_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    index_path = output_dir / "fasta_index.json"
    if index_path.exists():
        raise FileExistsError(f"Index already exists: {index_path}")

    with fasta_path.open("rb") as handle:
        header = handle.readline()
        if not header.startswith(b">"):
            raise ValueError("Expected a FASTA header beginning with >")
        sequence_offset = handle.tell()
        first_line = handle.readline()
        first_bases = first_line.rstrip(b"\r\n")
        if not first_bases:
            raise ValueError("FASTA contains no sequence bases")
        line_bases = len(first_bases)
        line_bytes = len(first_line)
        if line_bytes == line_bases:
            raise ValueError("Sequence lines must end in a newline")

        sequence_length = line_bases
        final_short_line = False
        for raw_line in handle:
            if raw_line.startswith(b">"):
                raise ValueError("Expected exactly one FASTA record")
            bases = raw_line.rstrip(b"\r\n")
            if not bases:
                if handle.read(1):
                    raise ValueError("Blank line inside the sequence")
                break
            if final_short_line:
                raise ValueError("Short sequence line before the final line")
            if len(bases) > line_bases:
                raise ValueError("Inconsistent FASTA line width")
            if len(bases) == line_bases and len(raw_line) != line_bytes:
                raise ValueError("Inconsistent FASTA newline width")
            if len(bases) < line_bases:
                final_short_line = True
            sequence_length += len(bases)

    if chunks > sequence_length:
        raise ValueError("More chunks than sequence bases")
    stat = fasta_path.stat()
    index = {
        "fasta_path": str(fasta_path),
        "description": header[1:].decode("ascii").strip(),
        "sequence_offset": sequence_offset,
        "line_bases": line_bases,
        "line_bytes": line_bytes,
        "sequence_length": sequence_length,
        "file_size": stat.st_size,
        "file_mtime_ns": stat.st_mtime_ns,
        "chunks": chunks,
        "context_bases": CONTEXT_BASES,
    }
    _write_json(index_path, index)
    return index


def load_index(output_dir: str | Path) -> dict:
    """Load the index and ensure the source FASTA has not changed."""
    with (Path(output_dir) / "fasta_index.json").open(encoding="utf-8") as handle:
        index = json.load(handle)
    stat = Path(index["fasta_path"]).stat()
    if stat.st_size != index["file_size"] or stat.st_mtime_ns != index["file_mtime_ns"]:
        raise ValueError("FASTA changed after it was indexed")
    return index


def core_bounds(sequence_length: int, chunks: int, chunk_id: int) -> tuple[int, int]:
    """Return the worker's zero-based, half-open core interval."""
    if not 0 <= chunk_id < chunks:
        raise ValueError("chunk_id outside indexed range")
    return (
        sequence_length * chunk_id // chunks,
        sequence_length * (chunk_id + 1) // chunks,
    )


def read_bases(index: dict, start: int, end: int) -> str:
    """Read a zero-based, half-open base range without FASTA line breaks."""
    if not 0 <= start <= end <= index["sequence_length"]:
        raise ValueError("Requested range is outside the sequence")
    if start == end:
        return ""
    width = index["line_bases"]
    stride = index["line_bytes"]
    sequence_offset = index["sequence_offset"]
    first_byte = sequence_offset + (start // width) * stride + start % width
    last_base = end - 1
    last_byte = (
        sequence_offset + (last_base // width) * stride + last_base % width
    )
    with Path(index["fasta_path"]).open("rb") as handle:
        handle.seek(first_byte)
        data = handle.read(last_byte + 1 - first_byte)
    sequence = data.replace(b"\n", b"").replace(b"\r", b"").decode("ascii")
    if len(sequence) != end - start:
        raise ValueError("FASTA data does not match the saved line index")
    return sequence.upper()


def extract_chunk(output_dir: str | Path, chunk_id: int) -> dict:
    """Write one core-plus-context sequence and its ownership metadata."""
    output_dir = Path(output_dir).resolve()
    index = load_index(output_dir)
    n = index["sequence_length"]
    core_start, core_end = core_bounds(n, index["chunks"], chunk_id)
    read_start = max(0, core_start - index["context_bases"])
    read_end = min(n, core_end + index["context_bases"])
    sequence = read_bases(index, read_start, read_end)

    stem = f"chunk_{chunk_id:02d}"
    sequence_path = output_dir / f"{stem}.seq"
    metadata_path = output_dir / f"{stem}.json"
    if sequence_path.exists() or metadata_path.exists():
        raise FileExistsError(f"Chunk {chunk_id} already exists in {output_dir}")
    metadata = {
        "chunk_id": chunk_id,
        "description": index["description"],
        "sequence_file": sequence_path.name,
        "core_start_1": core_start + 1,
        "core_end_1": core_end,
        "read_start_1": read_start + 1,
        "read_end_1": read_end,
        "core_offset_start": core_start - read_start,
        "core_offset_end_exclusive": core_end - read_start,
        "sequence_length": len(sequence),
    }

    temporary = sequence_path.with_name(f".{sequence_path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="ascii") as handle:
        handle.write(sequence)
    temporary.replace(sequence_path)
    _write_json(metadata_path, metadata)
    return metadata


def extract_local(output_dir: str | Path, workers: int) -> list[dict]:
    """Run independent chunk readers in local processes."""
    if workers < 1:
        raise ValueError("workers must be positive")
    index = load_index(output_dir)
    with ProcessPoolExecutor(max_workers=min(workers, index["chunks"])) as pool:
        return list(pool.map(
            extract_chunk,
            [str(output_dir)] * index["chunks"],
            range(index["chunks"]),
        ))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    index_command = commands.add_parser("index")
    index_command.add_argument("fasta")
    index_command.add_argument("output_dir")
    index_command.add_argument("--chunks", type=int, default=50)
    worker_command = commands.add_parser("worker")
    worker_command.add_argument("output_dir")
    worker_command.add_argument("chunk_id", type=int, nargs="?")
    local_command = commands.add_parser("local")
    local_command.add_argument("fasta")
    local_command.add_argument("output_dir")
    local_command.add_argument("--chunks", type=int, default=50)
    local_command.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()

    if args.command == "index":
        index = index_fasta(args.fasta, args.output_dir, args.chunks)
        print(
            f"Indexed {index['sequence_length']:,} bases "
            f"into {index['chunks']} chunks"
        )
    elif args.command == "worker":
        chunk_id = args.chunk_id
        if chunk_id is None:
            chunk_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
        metadata = extract_chunk(args.output_dir, chunk_id)
        print(
            f"Chunk {chunk_id}: owns {metadata['core_start_1']:,}–"
            f"{metadata['core_end_1']:,}; read {metadata['read_start_1']:,}–"
            f"{metadata['read_end_1']:,}"
        )
    else:
        index_fasta(args.fasta, args.output_dir, args.chunks)
        completed = extract_local(args.output_dir, args.workers)
        print(f"Wrote {len(completed)} sequence chunks to {args.output_dir}")


if __name__ == "__main__":
    main()
