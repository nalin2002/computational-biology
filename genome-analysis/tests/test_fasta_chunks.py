# AI_CODEX: Generated checks for parallel FASTA sequence extraction.
import json

from module_2_assignment.fasta_chunks import (
    core_bounds,
    extract_local,
    index_fasta,
    read_bases,
)


def test_parallel_chunks_preserve_bases_and_context(tmp_path):
    sequence = "ACGT" * 350 + "MRN"
    fasta_path = tmp_path / "chromosome.fasta"
    fasta_path.write_text(
        ">one chromosome\n"
        + "\n".join(sequence[i:i + 70] for i in range(0, len(sequence), 70))
        + "\n"
    )
    output_dir = tmp_path / "chunks"
    index = index_fasta(fasta_path, output_dir, chunks=4)
    assert index["sequence_length"] == len(sequence)
    assert read_bases(index, 69, 72) == sequence[69:72]

    completed = extract_local(output_dir, workers=2)
    assert len(completed) == 4
    reassembled = []
    for chunk_id in range(4):
        metadata = json.loads(
            (output_dir / f"chunk_{chunk_id:02d}.json").read_text()
        )
        stored = (output_dir / metadata["sequence_file"]).read_text()
        core_start, core_end = core_bounds(len(sequence), 4, chunk_id)
        read_start = max(0, core_start - 200)
        read_end = min(len(sequence), core_end + 200)
        assert stored == sequence[read_start:read_end]
        assert metadata["core_start_1"] == core_start + 1
        assert metadata["core_end_1"] == core_end
        reassembled.append(
            stored[metadata["core_offset_start"]:metadata["core_offset_end_exclusive"]]
        )
    assert "".join(reassembled) == sequence
