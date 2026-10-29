from module_2_assignment.dna_drive import decode_dna_to_text, encode_text_to_dna


def test_encode_decode_text_to_dna():

    # Test round-trip encoding and decoding of a simple string
    original_text = "Hello, World!"
    dna_sequence = encode_text_to_dna(original_text)
    decoded_text = decode_dna_to_text(dna_sequence)
    assert decoded_text == original_text

    # Test correct handling of uppercase, lowercase, digits, punctuation, and spaces
    test_string = "ABCDEFGHabcdu9 87!@#$"

    dna_sequence = encode_text_to_dna(test_string)
    decoded_string = decode_dna_to_text(dna_sequence)
    assert decoded_string == test_string

    # Test edge cases
    edge_cases = [ "","A","!" ]
    dna_sequences = [encode_text_to_dna(s) for s in edge_cases]
    decoded_strings = [decode_dna_to_text(dna) for dna in dna_sequences]
    assert decoded_strings == edge_cases