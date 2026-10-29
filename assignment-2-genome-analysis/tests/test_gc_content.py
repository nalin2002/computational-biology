from pytest import approx

from module_2_assignment.gc_content import gc_content


def test_gc_content():

    test_cases = {
        "ATCG": 50.0,
        "AAAA": 0.0,
        "CCCC": 100.0,
        "GGGG": 100.0,
        "TTTT": 0.0,
        "": 0.0,
        "ATGCAGTCAGTATGCGCGCGTACTGG": 57.69,
    }

    for sequence, expected in test_cases.items():
        assert gc_content(sequence) == approx(expected, abs=0.01 )


    
