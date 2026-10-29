import pytest
from module_3_assignment.hamming_dist import hamming_distance

def test_hamming_distance():

    # Test cases
    assert hamming_distance("1100111001", "0100010101") == 5
    assert hamming_distance("AGCTAAC", "AACTCCA") == 4
    assert hamming_distance("AGCTAAC", "AGCTAAC") == 0
    with pytest.raises(ValueError):
        hamming_distance("AGCTAAC", "AGCA")