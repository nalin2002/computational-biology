from module_2_assignment.cpg_islands import cpg_islands


def test_cpg_islands():
    test_cases = {
        "taacatacttattgtttttaactactcgttttccattcgactcatcacg"
        "ctccccccccccccccccccccttatccgttccgttcgacgtatttcgttgtct"
        "aatttctgacgtaacttgttccctgttaagtaccgtttatggcctatactccggt"
        "atttaaaacgacgacgattccaccgtaaagccgtcaaccagatgaacgacctcgc"
        "tcgttatatttttccggcaaaatccctatttccgattcgcttagtgctaccgacg"
        "ctatatcgttccgcaattcctcgagatcatcgatttcttctccggcgacgtctcaagt"
        "ttttccgttacaacgcgatctatcctgtaaattcgaccgcgctcattctcacgttt"
        "tatacattgcgcagttgattacgctaaataatccgctgactgttaccttccctgt"
        "tagattcgcgcattataaactacttactttaacaaacgattttcacagtttaat"
        "ttctgcgatgacgtctaactcttcagttttaaccgataacaaccttctcgacac"
        "ttcgtttcttataccatcctcgttatccatacccattcttaaatttctcac"
        "tactattctctttacaaccacattagctctaatcttacatctaatttctataca"
        "taaaatgctccttctgctgtatgtttctctttctcataattacatttttaatt"
        "actaaatccctcatccctcccacccatctattccaccatcaaggttatacac"
        "catgtattactgtaaaacccactaatattaattgtcaccgatattaaacgaa"
        "attcattcacacaaatttcattaattaccttttcttattaattgcatatgtac"
        "tctacatatactcaaccaactaaaaatcgatattttacatttgatttctaatg"
        "taccccacaactttcttgctttatgattgaacttagctttataataatagtta"
        "tttaccctaacgcatatactcttatccttatatgaaccttgcttatttgttag"
        "atttatccaatctaaaccacagataatatcccttctcttacttcattttatta"
        "tcaccattttcacttcttcctagatatatacaattatataactctattaccaca"
        "ttttcccttaactttctgttctgcactattatatttactctttttctaaaac"
        "cttcttaactttttcagatgca" : [(32,231),(33,232)]
    }

    for sequence, expected in test_cases.items():
        sequence = sequence.upper()
        results = cpg_islands(sequence)
        observed_islands= [start_end[:2] for start_end in results]

        print(f"Expected: {expected}, Observed: {observed_islands}")
        assert observed_islands == expected

