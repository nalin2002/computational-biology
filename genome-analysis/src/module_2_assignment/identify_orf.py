
class ProteinSyn:

    def __init__(self, description, sequence):
        self.description = description
        self.sequence = sequence.upper()

    # DNA -> RNA conversion- transcription
    def transcription(self, sequence):
        return sequence.replace("T", "U")

    # Generate reverse complement of DNA for generating the template 
    # strand and getting the reading frames
    def reverse_complement(self, sequence):

        watson_pairs = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "N": "N"
        }

        complement = ""

        for base in sequence:
            if base not in watson_pairs:
                continue
            complement += watson_pairs[base]

        return complement[::-1]


    # RNA -> Protein translation
    def translation(self, rna_sequence):

        codon_table = {
            "UUU": "F", "UUC": "F",
            "UUA": "L", "UUG": "L",
            "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
            "UAU": "Y", "UAC": "Y",
            "UAA": "*", "UAG": "*",
            "UGU": "C", "UGC": "C",
            "UGA": "*", "UGG": "W",

            "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
            "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAU": "H", "CAC": "H",
            "CAA": "Q", "CAG": "Q",
            "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",

            "AUU": "I", "AUC": "I", "AUA": "I",
            "AUG": "M",
            "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "AAU": "N", "AAC": "N",
            "AAA": "K", "AAG": "K",
            "AGU": "S", "AGC": "S",
            "AGA": "R", "AGG": "R",

            "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
            "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GAU": "D", "GAC": "D",
            "GAA": "E", "GAG": "E",
            "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
        }

        protein = ""

        for i in range(0, len(rna_sequence) - 2, 3):

            # getting non overlapping groups of 3 letters
            codon = rna_sequence[i:i + 3]

            # not including stop codon in protein
            if codon in ("UAA", "UAG", "UGA"):
                break

            if codon not in codon_table:
                continue

            protein += codon_table[codon]

        return protein


    # Finding open reading frames(ORF) in one DNA strand
    def find_orfs_in_strand(self, dna_sequence, frame, reverse=False):

        stop_codons = {"TAA", "TAG", "TGA"}
        orfs = []

        # 3 frames with 3 offsets
        offset = abs(frame) - 1

        for i in range(offset, len(dna_sequence) - 2, 3):
            codon = dna_sequence[i:i + 3]

            # Every ATG begins a potential ORF
            if codon == "ATG":

                # Search for first stop codon in same frame
                for j in range(i + 3, len(dna_sequence) - 2, 3):
                    stop = dna_sequence[j:j + 3]
                    if stop in stop_codons:

                        # DNA from ATG up to the stop codon, excluding the stop codon
                        coding_dna = dna_sequence[i:j]

                        # DNA -> RNA
                        rna = self.transcription(coding_dna)

                        # RNA -> protein
                        protein = self.translation(rna)

                        if not reverse:
                            # indexing with 1-based coordinates
                            start_position = i + 1

                            # j is first base of stop codon, and +3 gives its
                            #  final base in 1-based indexing.
                            stop_position = j + 3

                        else:
                            # Convert reverse-complement positions
                            # back to coordinates on original sequence
                            n = len(self.sequence)

                            start_position = n - i
                            stop_position = n - j - 2

                        orfs.append({
                            "frame": frame,
                            "start": start_position,
                            "stop": stop_position,
                            "length": len(protein),
                            "protein": protein
                        })

                        # ORF ends at first break
                        break

        return orfs


    def find_all_orfs(self):

        all_orfs = []

        # in forward strands
        for frame in [1, 2, 3]:

            orfs = self.find_orfs_in_strand(
                self.sequence,
                frame,
                reverse=False
            )

            all_orfs.extend(orfs)

        # getting the template strand
        reverse_seq = self.reverse_complement(self.sequence)

        # going from the back
        for frame in [-1, -2, -3]:

            orfs = self.find_orfs_in_strand(
                reverse_seq,
                frame,
                reverse=True
            )

            all_orfs.extend(orfs)

        return all_orfs

    def find_orfs_postive_strand(self):

        all_orfs=[]

        for frame in [1, 2, 3]:

            orfs = self.find_orfs_in_strand(
                self.sequence,
                frame,
                reverse=False
            )

            all_orfs.extend(orfs)

        return all_orfs

    def find_orfs_negative_strand(self):

        all_orfs=[]

        reverse_seq = self.reverse_complement(self.sequence)

        for frame in [-1, -2, -3]:

            orfs = self.find_orfs_in_strand(
                reverse_seq,
                frame,
                reverse=True
            )

            all_orfs.extend(orfs)

        return all_orfs