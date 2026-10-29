from module_2_assignment.gene_music import base_mapping_frequencies


def test_music_output():
    # Testing nucleotide-to-note mapping returns the correct frequency for each base
    assert base_mapping_frequencies["A"] == 523.25
    assert base_mapping_frequencies["T"] == 587.33
    assert base_mapping_frequencies["G"] == 659.25
    assert base_mapping_frequencies["C"] == 783.99
    