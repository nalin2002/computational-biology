# AI_CODEX: Generated these integration checks for the parallel pipeline.
"""Integration checks for chromosome chunk boundaries and the FASTA index."""

import json

from module_2_assignment.cpg_islands import cpg_islands
from module_2_assignment.gene_finder import merge_cpg_islands
from module_2_assignment.identify_orf import ProteinSyn
from module_2_assignment.parallel_gene_finder import (
    filter_chunk,
    finalize,
    index_fasta,
    prepare,
    read_bases,
    scan_chunk,
)


def reference_candidates(sequence):
    """Use the original ORF/CpG functions, with the requested CpG merge."""
    islands = merge_cpg_islands(cpg_islands(sequence))
    orfs = ProteinSyn("test", sequence)
    results = []
    for strand, records in ((1, orfs.find_orfs_postive_strand()),
                            (-1, orfs.find_orfs_negative_strand())):
        for orf in records:
            start, stop = orf["start"], orf["stop"]
            length = abs(stop - start) + 1
            if length < 96:
                continue
            if strand == 1:
                first_lo, first_hi = start, start + length // 3
            else:
                first_lo, first_hi = start - length // 3, start
            if any(a <= first_hi and b >= first_lo for a, b in islands):
                results.append([start, stop, length, strand])
    return results


def test_parallel_pipeline_crosses_both_orf_and_cpg_boundaries(tmp_path):
    plus = "CG" * 20 + "ATG" + "CGC" * 31 + "TAA" + "A" * 100
    reverse = plus.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    sequence = plus + "N" * 20 + reverse
    fasta = tmp_path / "example.fasta"
    with fasta.open("w") as handle:
        handle.write(">chromosome\n")
        for i in range(0, len(sequence), 70):
            handle.write(sequence[i:i + 70] + "\n")
    output_dir = tmp_path / "out"
    index_fasta(fasta, output_dir, tasks=7)
    meta = json.loads((output_dir / "index.json").read_text())
    assert read_bases(meta, 68, 73) == sequence[68:73]
    assert read_bases(meta, 0, len(sequence)) == sequence

    for task_id in range(7):
        scan_chunk(output_dir, task_id)
    prepare(output_dir)
    for task_id in range(7):
        filter_chunk(output_dir, task_id)
    finalize(output_dir)

    rows = (output_dir / "gene_candidates.tsv").read_text().splitlines()[1:]
    observed = [[int(start), int(stop), int(length), 1 if sign == "+" else -1]
                for start, stop, length, sign in (row.split("\t") for row in rows)]
    expected = reference_candidates(sequence)
    assert observed == expected
    assert [41, 139, 99, 1] in observed
    assert [458, 360, 99, -1] in observed
    summary = json.loads((output_dir / "summary.json").read_text())
    assert summary["gene_candidates"] == len(expected)
    assert summary["merged_cpg_islands"] == len(
        merge_cpg_islands(cpg_islands(sequence))
    )


def test_ambiguous_base_keeps_reverse_coordinates(tmp_path):
    plus = "CG" * 20 + "ATG" + "CGC" * 31 + "TAA" + "A" * 100
    reverse = plus.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    sequence = plus + "N" * 10 + "M" + "N" * 10 + reverse
    fasta = tmp_path / "ambiguous.fasta"
    fasta.write_text(">chromosome\n" + "\n".join(
        sequence[i:i + 70] for i in range(0, len(sequence), 70)
    ) + "\n")
    output_dir = tmp_path / "out"
    index_fasta(fasta, output_dir, tasks=5)
    for task_id in range(5):
        scan_chunk(output_dir, task_id)
    prepare(output_dir)
    for task_id in range(5):
        filter_chunk(output_dir, task_id)
    finalize(output_dir)
    rows = (output_dir / "gene_candidates.tsv").read_text()
    assert "459\t361\t99\t-\n" in rows
