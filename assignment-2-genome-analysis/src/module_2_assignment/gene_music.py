import numpy as np
from scipy.io.wavfile import write

# Frequencies for individual bases
base_mapping_frequencies = {
    'A': 523.25,   # C5
    'T': 587.33,   # D5
    'G': 659.25,   # E5
    'C': 783.99    # G5
}

# Frequencies for amino acids based on codon usage
amino_acid_codon_mapping_frequencies = {
    'Met': 65.41,    # C2
    'Gly': 73.42,    # D2
    'Ala': 82.41,    # E2
    'Val': 98.00,    # G2
    'Leu': 110.00,   # A2

    'Ile': 123.47,   # B2
    'Pro': 130.81,   # C3
    'Phe': 146.83,   # D3
    'Tyr': 164.81,   # E3
    'Trp': 196.00,   # G3

    'Ser': 220.00,   # A3
    'Thr': 246.94,   # B3
    'Cys': 261.63,   # C4
    'Asn': 293.66,   # D4
    'Gln': 329.63,   # E4

    'Lys': 392.00,   # G4
    'Arg': 440.00,   # A4
    'His': 493.88,   # B4
    'Asp': 523.25,   # C5
    'Glu': 587.33,   # D5

    'Stop': 65.41    # C2
}

codon_table = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop','UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop','UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

# AI_Claude: Used Claude just to help improve the code for 
# creatively generating music from DNA sequences.
def generate_music(sequence):
    sequence = sequence.upper()[:300]

    base_notes = [
    base_mapping_frequencies[base]
    for base in sequence
    if base in base_mapping_frequencies]

    amino_acid_notes = {}

    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]

        if len(codon) != 3:
            break

        amino_acid = codon_table.get(codon, 'Stop')

        freq = amino_acid_codon_mapping_frequencies.get(
            amino_acid, 65.41
        )

        amino_acid_notes[i // 3] = freq


    sample_rate = 44100

    # Faster, more rhythmic
    base_duration = 0.12

    num_samples = int(sample_rate * base_duration)

    audio = np.zeros(
        num_samples * len(base_notes)
    )


    for i, freq in enumerate(base_notes):

        t = np.linspace(
            0,
            base_duration,
            num_samples,
            False
        )

        # Very fast attack + exponential decay
        pluck_envelope = np.exp(-25 * t)

        # Fundamental + harmonics gives a brighter sound
        base_tone = (
            0.30 * np.sin(2 * np.pi * freq * t) +
            0.12 * np.sin(2 * np.pi * freq * 2 * t) +
            0.05 * np.sin(2 * np.pi * freq * 3 * t)
        )

        base_tone *= pluck_envelope

        codon_idx = i // 3

        aa_tone = np.zeros_like(t)

        # Trigger AA percussion only on first base of codon
        if i % 3 == 0 and codon_idx in amino_acid_notes:

            aa_freq = amino_acid_notes[codon_idx]

            # Much faster decay → drum-like
            percussion_envelope = np.exp(-35 * t)

            aa_tone = (
                0.35 * np.sin(2 * np.pi * aa_freq * t) +
                0.15 * np.sin(2 * np.pi * aa_freq * 2 * t)
            )

            aa_tone *= percussion_envelope

        kick = np.zeros_like(t)

        # Kick every 6 bases = every 2 codons
        if i % 6 == 0:

            kick_envelope = np.exp(-30 * t)

            # Falling-frequency kick
            kick_freq = 60 + 60 * np.exp(-40 * t)

            phase = 2 * np.pi * np.cumsum(kick_freq) / sample_rate

            kick = (
                0.45 *
                np.sin(phase) *
                kick_envelope
            )

        clap = np.zeros_like(t)

        # Between kicks
        if i % 6 == 3:

            noise = np.random.normal(0, 1, len(t))

            clap_envelope = np.exp(-45 * t)

            clap = (
                0.10 *
                noise *
                clap_envelope
            )

        tone = (
            base_tone +
            aa_tone +
            kick +
            clap
        )

        start_idx = i * num_samples
        end_idx = start_idx + num_samples

        audio[start_idx:end_idx] += tone

    max_val = np.max(np.abs(audio))

    if max_val > 0:
        audio = audio / max_val


    write('../data/gene_bdnf.wav', sample_rate, (audio * 32767).astype(np.int16) )

    print("DNA music generated successfully!")

    return base_notes, amino_acid_notes, audio