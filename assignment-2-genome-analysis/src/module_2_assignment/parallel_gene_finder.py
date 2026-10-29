# AI_CODEX: Generated this parallel pipeline and its boundary handling.
"""Parallel chromosome gene finder with exact cross-chunk ORF and CpG handling.

Run the stages through slurm/submit_gene_array.sh or directly with
``uv run python -m module_2_assignment.parallel_gene_finder <stage> ...``.
All chromosome coordinates in intermediate files are zero-based; output is
one-based and follows the existing gene_finder.py conventions.
"""

import argparse
import bisect
import json
import os
import time
from pathlib import Path

from module_2_assignment.cpg_islands import cpg_islands
from module_2_assignment.gene_finder import merge_cpg_islands

STOP_CODONS = {"TAA", "TAG", "TGA"}
# Ambiguous symbols cannot form an exact start or stop codon. Treat them as N
# on the reverse strand so chromosome coordinates remain unchanged.
COMPLEMENT = str.maketrans("ACGTNMR", "TGCANNN")


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    with temporary.open("w") as handle:
        json.dump(value, handle, separators=(",", ":"))
    temporary.replace(path)


def read_json(path):
    with Path(path).open() as handle:
        return json.load(handle)


def index_fasta(fasta_path, output_dir, tasks):
    """Validate and index one uncompressed, fixed-width FASTA record."""
    started = time.time()
    fasta_path = Path(fasta_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if (output_dir / "index.json").exists():
        raise FileExistsError(f"Output directory already contains a run: {output_dir}")
    if tasks < 1 or tasks > 50:
        raise ValueError("tasks must be between 1 and 50")

    with fasta_path.open("rb") as handle:
        header = handle.readline()
        if not header.startswith(b">"):
            raise ValueError("FASTA must begin with a > header")
        sequence_offset = handle.tell()
        first = handle.readline()
        first_bases = first.rstrip(b"\r\n")
        if not first_bases:
            raise ValueError("FASTA sequence is empty")
        line_bases = len(first_bases)
        line_bytes = len(first)
        if line_bytes == line_bases:
            raise ValueError("Sequence lines must have newlines")
        length = 0
        short_line_seen = False
        length += len(first_bases)
        for raw_line in handle:
            if raw_line.startswith(b">"):
                raise ValueError("This workflow expects exactly one FASTA record")
            line = raw_line.rstrip(b"\r\n")
            if not line:
                # A single empty line after the last sequence line is accepted.
                if handle.read(1):
                    raise ValueError("Blank line inside the sequence")
                break
            if short_line_seen:
                raise ValueError("Short sequence line before the final line")
            if len(line) > line_bases or len(raw_line) > line_bytes:
                raise ValueError("Inconsistent FASTA line width")
            if len(line) == line_bases and len(raw_line) != line_bytes:
                raise ValueError("Inconsistent FASTA newline width")
            if len(line) < line_bases:
                short_line_seen = True
            length += len(line)

    size = fasta_path.stat().st_size
    write_json(output_dir / "index.json", {
        "fasta": str(fasta_path),
        "description": header[1:].decode("ascii").strip(),
        "sequence_offset": sequence_offset,
        "line_bases": line_bases,
        "line_bytes": line_bytes,
        "length": length,
        "size": size,
        "mtime_ns": fasta_path.stat().st_mtime_ns,
        "tasks": tasks,
        "started_at": started,
    })
    print(f"Indexed {length:,} bases in {tasks} chunks")


def load_index(output_dir):
    meta = read_json(Path(output_dir) / "index.json")
    stat = Path(meta["fasta"]).stat()
    if stat.st_size != meta["size"] or stat.st_mtime_ns != meta["mtime_ns"]:
        raise ValueError("FASTA changed after indexing")
    return meta


def chunk_bounds(length, tasks, task_id):
    if not 0 <= task_id < tasks:
        raise ValueError("task ID outside indexed range")
    return length * task_id // tasks, length * (task_id + 1) // tasks


def read_bases(meta, start, end):
    """Read a zero-based half-open genomic range, ignoring FASTA newlines."""
    if not 0 <= start <= end <= meta["length"]:
        raise ValueError("Base interval outside FASTA sequence")
    if start == end:
        return ""
    width = meta["line_bases"]
    stride = meta["line_bytes"]
    offset = meta["sequence_offset"]
    first_byte = offset + (start // width) * stride + start % width
    last_base = end - 1
    last_byte = offset + (last_base // width) * stride + last_base % width
    with Path(meta["fasta"]).open("rb") as handle:
        handle.seek(first_byte)
        data = handle.read(last_byte + 1 - first_byte)
    sequence = data.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
    if len(sequence) != end - start:
        raise ValueError("FASTA index does not match sequence layout")
    return sequence


def scan_chunk(output_dir, task_id):
    """Own codons and CpG windows by their leftmost genomic base."""
    started = time.perf_counter()
    meta = load_index(output_dir)
    lo, hi = chunk_bounds(meta["length"], meta["tasks"], task_id)
    # Two extra bases complete codons; 199 complete 200-base CpG windows.
    sequence = read_bases(meta, lo, min(meta["length"], hi + 199))
    starts = [[[] for _ in range(3)] for _ in range(2)]
    stops = [[[] for _ in range(3)] for _ in range(2)]
    for p in range(lo, min(hi, meta["length"] - 2)):
        codon = sequence[p - lo:p - lo + 3]
        frame = p % 3
        if codon == "ATG":
            starts[0][frame].append(p)
        elif codon in STOP_CODONS:
            stops[0][frame].append(p)
        reverse_codon = codon.translate(COMPLEMENT)[::-1]
        if reverse_codon == "ATG":
            starts[1][frame].append(p)
        elif reverse_codon in STOP_CODONS:
            stops[1][frame].append(p)

    windows = cpg_islands(sequence)
    # cpg_islands returns one-based coordinates within its input. A window
    # starting in the 199-base context belongs to the next chunk.
    owned_windows = [(lo + start, lo + end) for start, end, _, _ in windows
                     if lo + start - 1 < hi]
    islands = merge_cpg_islands(owned_windows)
    write_json(Path(output_dir) / f"scan_{task_id:02d}.json", {
        "task": task_id, "lo": lo, "hi": hi,
        "starts": starts, "stops": stops, "islands": islands,
        "seconds": time.perf_counter() - started,
    })
    print(f"Chunk {task_id}: {lo + 1:,}-{hi:,}; {len(islands):,} local CpG regions")


def prepare(output_dir):
    """Merge CpG regions globally and build six frame-specific stop indexes."""
    started = time.perf_counter()
    meta = load_index(output_dir)
    all_islands = []
    all_stops = [[[] for _ in range(3)] for _ in range(2)]
    scan_seconds = []
    for task_id in range(meta["tasks"]):
        chunk = read_json(Path(output_dir) / f"scan_{task_id:02d}.json")
        lo, hi = chunk_bounds(meta["length"], meta["tasks"], task_id)
        if (chunk["task"], chunk["lo"], chunk["hi"]) != (task_id, lo, hi):
            raise ValueError(f"Chunk {task_id} has incorrect ownership")
        all_islands.extend(chunk["islands"])
        scan_seconds.append(chunk["seconds"])
        for strand in range(2):
            for frame in range(3):
                all_stops[strand][frame].extend(chunk["stops"][strand][frame])
    merged = merge_cpg_islands(all_islands)
    write_json(Path(output_dir) / "prepared.json", {
        "islands": merged, "stops": all_stops,
        "scan_max_seconds": max(scan_seconds),
        "scan_sum_seconds": sum(scan_seconds),
        "prepare_seconds": time.perf_counter() - started,
    })
    print(f"Global CpG regions after merge: {len(merged):,}")


def filter_chunk(output_dir, task_id):
    """Pair each owned ATG with the first downstream in-frame stop."""
    started = time.perf_counter()
    meta = load_index(output_dir)
    if not 0 <= task_id < meta["tasks"]:
        raise ValueError("task ID outside indexed range")
    prepared = read_json(Path(output_dir) / "prepared.json")
    chunk = read_json(Path(output_dir) / f"scan_{task_id:02d}.json")
    islands = prepared["islands"]
    island_ends = [end for _, end in islands]
    candidates = []
    for strand in range(2):
        for frame in range(3):
            stop_positions = prepared["stops"][strand][frame]
            for p in chunk["starts"][strand][frame]:
                if strand == 0:
                    i = bisect.bisect_right(stop_positions, p)
                    if i == len(stop_positions):
                        continue
                    q = stop_positions[i]
                    start, stop = p + 1, q + 3
                    length = stop - start + 1
                    first_lo, first_hi = start, start + length // 3
                    sign = 1
                else:
                    i = bisect.bisect_left(stop_positions, p) - 1
                    if i < 0:
                        continue
                    q = stop_positions[i]
                    start, stop = p + 3, q + 1
                    length = start - stop + 1
                    first_lo, first_hi = start - length // 3, start
                    sign = -1
                if length < 96:
                    continue
                j = bisect.bisect_left(island_ends, first_lo)
                if j < len(islands) and islands[j][0] <= first_hi:
                    candidates.append([start, stop, length, sign])
    write_json(Path(output_dir) / f"candidates_{task_id:02d}.json", {
        "task": task_id, "candidates": candidates,
        "seconds": time.perf_counter() - started,
    })
    print(f"Chunk {task_id}: {len(candidates):,} gene candidates")


def finalize(output_dir):
    started = time.perf_counter()
    meta = load_index(output_dir)
    prepared = read_json(Path(output_dir) / "prepared.json")
    candidates = []
    filter_seconds = []
    for task_id in range(meta["tasks"]):
        part = read_json(Path(output_dir) / f"candidates_{task_id:02d}.json")
        if part["task"] != task_id:
            raise ValueError(f"Candidate chunk {task_id} has wrong task ID")
        candidates.extend(part["candidates"])
        filter_seconds.append(part["seconds"])
    # The sequential implementation reports all positive candidates first,
    # then negative candidates in reverse genomic order.
    candidates.sort(key=lambda row: (0, row[0]) if row[3] == 1
                    else (1, -row[0]))
    output_dir = Path(output_dir)
    with (output_dir / "gene_candidates.tsv").open("w") as handle:
        handle.write("start\tstop\tlength_bp\tstrand\n")
        for start, stop, length, strand in candidates:
            handle.write(f"{start}\t{stop}\t{length}\t{'+' if strand == 1 else '-'}\n")
    with (output_dir / "gene_candidates.txt").open("w") as handle:
        for start, stop, length, strand in candidates:
            handle.write(f"Gene candidate: position {start} to {stop} "
                         f"(strand: {'+' if strand == 1 else '-'}, "
                         f"length: {length} bp)\n")
        handle.write(f"Total number of gene candidates found: {len(candidates)}\n")
    summary = {
        "description": meta["description"],
        "sequence_bases": meta["length"],
        "tasks": meta["tasks"],
        "merged_cpg_islands": len(prepared["islands"]),
        "gene_candidates": len(candidates),
        "pipeline_wall_seconds": time.time() - meta["started_at"],
        "scan_max_seconds": prepared["scan_max_seconds"],
        "scan_sum_seconds": prepared["scan_sum_seconds"],
        "prepare_seconds": prepared["prepare_seconds"],
        "filter_max_seconds": max(filter_seconds),
        "filter_sum_seconds": sum(filter_seconds),
        "finalize_seconds": time.perf_counter() - started,
    }
    write_json(output_dir / "summary.json", summary)
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="stage", required=True)
    index = sub.add_parser("index")
    index.add_argument("fasta")
    index.add_argument("output_dir")
    index.add_argument("--tasks", type=int, default=50)
    for name in ("scan", "filter"):
        command = sub.add_parser(name)
        command.add_argument("output_dir")
        command.add_argument("task_id", type=int, nargs="?")
    for name in ("prepare", "finalize"):
        sub.add_parser(name).add_argument("output_dir")
    args = parser.parse_args()
    if args.stage == "index":
        index_fasta(args.fasta, args.output_dir, args.tasks)
    elif args.stage == "scan":
        task_id = (args.task_id if args.task_id is not None
                   else int(os.environ["SLURM_PROCID"]))
        scan_chunk(args.output_dir, task_id)
    elif args.stage == "prepare":
        prepare(args.output_dir)
    elif args.stage == "filter":
        task_id = (args.task_id if args.task_id is not None
                   else int(os.environ["SLURM_PROCID"]))
        filter_chunk(args.output_dir, task_id)
    else:
        finalize(args.output_dir)


if __name__ == "__main__":
    main()
