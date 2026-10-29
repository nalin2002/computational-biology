from module_3_assignment.matrices import PAM250, BLOSUM62
import argparse
from module_3_assignment.fasta_parser import read_fasta
import numpy as np

# The application should read sequences from FASTA-formatted files and allow the user to select the type of alignment (local or global), the substitution matrix, and the gap penalty from the command line.

# **Supported substitution matrices:** BLOSUM62 and PAM250. Python dictionaries for both matrices are located in the `src/` module.

# **Gap penalty:** A linear (per-gap-character) penalty applied as a negative value. Pass it as a negative number on the command line (e.g., `--gap_penalty -2` means each gap character costs −2).

# An example of running the program:

# ```bash
# python align.py --seq1 data/sequence1.fasta \
#                 --seq2 data/sequence2.fasta \
#                 --type global \
#                 --matrix blosum62 \
#                 --gap_penalty -2
# ```

# **Output format:** The program should report the alignment, the alignment score, and the sequence identity. Sequence identity is defined as the number of matched (identical) positions divided by the total alignment length (including gap columns), expressed as a percentage.

# As a sample test case, align human **Bradykinin** (`RPPGFSPFR`) against human **Kallidin** (`KRPPGFSPFR`) — two related vasoactive peptides that differ by a single N-terminal lysine residue. Using BLOSUM62 and a gap penalty of −2, the expected output is:

# ```
# seq 1   -RPPGFSPFR 
# seq 2   KRPPGFSPFR 
# Alignment Score:    51
# Sequence Identity:  90%
# ```

# The FASTA files for these sequences are included in `data/` as `bradykinin.fasta` and `kallidin.fasta`.

# Your program should exit gracefully with a helpful error message if given invalid input (e.g., non-FASTA files, unrecognised matrix name, non-protein characters in sequences).

def global_alignment(seq1, seq2, matrix, gap_penalty):
    n = len(seq1)
    m = len(seq2)
    # Initialize the scoring matrix
    score_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # Filling the first row and column with gap penalties
    for i in range(1, n + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] + gap_penalty
    for j in range(1, m + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] + gap_penalty

    # Fill the scoring matrix
    for i in range(1,n+1):
        for j in range(1,m+1):
            gap_j = score_matrix[i][j - 1] + gap_penalty
            gap_i = score_matrix[i - 1][j] + gap_penalty
            match= score_matrix[i - 1][j - 1] + matrix[seq1[i - 1]][seq2[j - 1]]
            score_matrix[i][j] = max(gap_i, gap_j, match)

    alignment_score = score_matrix[n][m]

    # Backtrack to find the aligned sequences
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = n, m
    while i > 0 or j > 0:
        current_score = score_matrix[i][j]
        if i > 0 and j > 0 and current_score == score_matrix[i - 1][j - 1] + matrix[seq1[i - 1]][seq2[j - 1]]:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and current_score == score_matrix[i - 1][j] + gap_penalty:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2, alignment_score

def local_alignment(seq1, seq2, matrix, gap_penalty):
    n = len(seq1)
    m = len(seq2)
    # Initialize the scoring matrix
    score_matrix = np.zeros((n + 1, m + 1), dtype=int)
    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + matrix[seq1[i - 1]][seq2[j - 1]]
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(0,match, delete, insert)

    alignment_score = np.max(score_matrix)

    # Backtrack to find the aligned sequences
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)
    while i > 0 and j > 0 and score_matrix[i][j] != 0:
        current_score = score_matrix[i][j]
        if current_score == score_matrix[i - 1][j - 1] + matrix[seq1[i - 1]][seq2[j - 1]]:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif current_score == score_matrix[i - 1][j] + gap_penalty:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1
    return aligned_seq1, aligned_seq2, alignment_score

def main():
    # to execute this file from the command line inside , you can use the following command: 
    # first cd src/module_3_assignment   
    # uv run align.py --seq1 ../../data/bradykinin.fasta --seq2 ../../data/kallidin.fasta --type global --matrix blosum62 --gap_penalty -2
    parser=argparse.ArgumentParser(description="Align two protein sequences.")
    
    parser.add_argument("--seq1", required=True, help="Path to the first FASTA file.")
    parser.add_argument("--seq2", required=True, help="Path to the second FASTA file.")
    parser.add_argument("--type", required=True, choices=["global", "local"], help="Type of alignment: 'global' or 'local'.")
    parser.add_argument("--matrix", required=True, choices=["blosum62", "pam250"], help="Substitution matrix: 'blosum62' or 'pam250'.")
    parser.add_argument("--gap_penalty", required=True, type=int, help="Gap penalty (negative integer).")

    args = parser.parse_args()
    sequence1 =""
    sequence2 =""
    # Read sequences from FASTA files
    for description, sequence in read_fasta(args.seq1):
        seq1_description = description
        sequence1 = sequence

    for description, sequence in read_fasta(args.seq2):
        seq2_description = description
        sequence2 = sequence

    # Select substitution matrix
    if args.matrix == "blosum62":
        matrix = BLOSUM62
    elif args.matrix == "pam250":
        matrix = PAM250
    else:
        raise ValueError("Invalid substitution matrix. Choose 'blosum62' or 'pam250'.")

    # Perform alignment
    if args.type == "global":  
        aligned_seq1, aligned_seq2, score  = global_alignment(sequence1, sequence2, matrix, args.gap_penalty)

        seq_identity=0
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == b:
                seq_identity += 1
        seq_identity_percentage = (seq_identity / len(aligned_seq1)) * 100

        print(f"seq 1   {aligned_seq1}")
        print(f"seq 2   {aligned_seq2}")
        print(f"Alignment Score:    {score}")
        print(f"Sequence Identity:  {seq_identity_percentage:.0f}%")

    elif args.type == "local":
        aligned_seq1, aligned_seq2, score  = local_alignment(sequence1, sequence2, matrix, args.gap_penalty)

        seq_identity=0
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == b:
                seq_identity += 1
        seq_identity_percentage = (seq_identity / len(aligned_seq1)) * 100

        print(f"seq 1   {aligned_seq1}")
        print(f"seq 2   {aligned_seq2}")
        print(f"Alignment Score:    {score}")
        print(f"Sequence Identity:  {seq_identity_percentage:.0f}%")
    else:
        raise ValueError("Invalid alignment type. Choose 'global' or 'local'.")

    


if __name__ == "__main__":
    main()