from module_2_assignment.gene_finder import find_genes


def test_find_genes():

    # Test case 1: ORF >= 96 bp with a CpG island overlapping the first third
    dna_sequence_1 = (
        "CG" * 20
        + "ATG"
        + "CGC" * 31
        + "TAA"
        + "A" * 100
    )
    description_1 = "Test sequence 1"
    expected_output_1 = [ [41, 139, 99, 1] ]
    assert find_genes(dna_sequence_1, description_1) == expected_output_1

    # Test case 2: Short ORF (< 96 bp) even though it overlaps a CpG-rich region
    dna_sequence_2 = (
        "CG" * 20
        + "ATG"
        + "CGC" * 10
        + "TAA"
        + "CG" * 70
    )
    description_2 = "Test sequence 2"
    expected_output_2 = []
    assert find_genes(dna_sequence_2, description_2) == expected_output_2


    # Test case 3: ORF >= 96 bp but no overlapping CpG island
    dna_sequence_3 = (
        "A" * 40
        + "ATG"
        + "AAA" * 31
        + "TAA"
        + "A" * 100
    )
    description_3 = "Test sequence 3"
    expected_output_3 = []
    assert find_genes(dna_sequence_3, description_3) == expected_output_3