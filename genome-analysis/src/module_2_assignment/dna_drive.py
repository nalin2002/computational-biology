
def encode_text_to_dna(text):

    # Convert text to binary representation
    binary_string = ''.join(format(ord(char), '08b') for char in text)

    # Map binary to DNA bases
    dna_mapping = {'00': 'A', '01': 'T', '10': 'G', '11': 'C'}
    dna_sequence = ''.join(dna_mapping[binary_string[i:i + 2]] 
                           for i in range(0, len(binary_string), 2))

    # each character is represented by 8 bits, which corresponds to 4 DNA bases
    # So, the length of the DNA sequence should be 4 times the length of the input text

    return dna_sequence


def decode_dna_to_text(dna_sequence):

    # Map DNA bases back to binary
    binary_mapping = {'A': '00', 'T': '01', 'G': '10', 'C': '11'}
    binary_string = ''.join(binary_mapping[base] for base in dna_sequence)

    # Convert binary representation back to text
    text = ''.join(chr(int(binary_string[i:i + 8], 2)) 
                   for i in range(0, len(binary_string), 8))

    return text


