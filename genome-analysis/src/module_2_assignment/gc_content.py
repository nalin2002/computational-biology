
def gc_content(sequence: str) -> float:
    """
    Calculate the GC-content of a DNA sequence.
    Args:
        sequence: A string of nucleotide characters (A, T, G, C).
    Returns:
        The GC-content as a percentage (e.g., 50.0 for a sequence
        that is half G/C).
    """

    sequence = sequence.upper()

    # we just need to calculate the number of G and C characters in the sequence

    gc_count = sum(1 for char in sequence if char in "GC")

    # If the sequence is empty, we return 0.0 to avoid division by zero
    if len(sequence) == 0:
        return 0.0

    # Calculating the GC-content as a percentage
    return (gc_count / len(sequence)) * 100.0

