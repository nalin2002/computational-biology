from module_2_assignment.cpg_islands import cpg_islands
from module_2_assignment.identify_orf import ProteinSyn


def merge_cpg_islands(cpg_islands_list):

    if not cpg_islands_list:
        return []

    merged = []

    current_start = cpg_islands_list[0][0]
    current_end = cpg_islands_list[0][1]

    for cpg_island in cpg_islands_list[1:]:

        cpg_start = cpg_island[0]
        cpg_end = cpg_island[1]

        # Current CpG window overlaps the previous region
        if cpg_start <= current_end:
            current_end = max(current_end, cpg_end)

        else:
            merged.append((current_start, current_end))

            # Starting a new region
            current_start = cpg_start
            current_end = cpg_end

    # Adding the final region
    merged.append((current_start, current_end))

    return merged

def find_genes(dna_sequence, description):

    print("Starting gene finder...", flush=True)

    protein_syn = ProteinSyn(description, dna_sequence)

    # Finding ORFs from positive strand
    positive_orfs = protein_syn.find_orfs_postive_strand()
    print(f"Positive strand ORFs found: {len(positive_orfs):,}", flush=True)

    # Finding ORFs from negative strand
    negative_orfs = protein_syn.find_orfs_negative_strand()
    print(f"Negative strand ORFs found: {len(negative_orfs):,}", flush=True)

    # Identifying CpG islands
    cpg_islands_list = cpg_islands(dna_sequence)
    cpg_islands_list = merge_cpg_islands(cpg_islands_list)
    print(f"CpG islands found: {len(cpg_islands_list):,}", flush=True)

    genes_candidates = []

    # Positive strand
    for orf in positive_orfs:

        orf_start = orf["start"]
        orf_end = orf["stop"]

        orf_length_bp = abs(orf_end - orf_start) + 1

        if orf_length_bp >= 96:

            # First third of ORF
            first_third_end = orf_start + (orf_length_bp // 3)

            for cpg_island in cpg_islands_list:

                cpg_start = cpg_island[0]
                cpg_end = cpg_island[1]

                # Check if CpG island overlaps first third of ORF
                if (cpg_start <= first_third_end) and (cpg_end >= orf_start):

                    genes_candidates.append(
                        [orf_start, orf_end, orf_length_bp, 1]
                    )
                    break


    # Negative strand
    for orf in negative_orfs:

        # Negative strand start is the higher genomic coordinate
        orf_start = orf["start"]
        orf_end = orf["stop"]

        orf_length_bp = orf_start - orf_end + 1

        if orf_length_bp >= 96:

            # Negative strand proceeds toward lower coordinates
            first_third_end = orf_start - (orf_length_bp // 3)

            for cpg_island in cpg_islands_list:

                cpg_start = cpg_island[0]
                cpg_end = cpg_island[1]

                # First third is [first_third_end, orf_start]
                if (cpg_start <= orf_start) and (cpg_end >= first_third_end):

                    genes_candidates.append(
                        [orf_start, orf_end, orf_length_bp, -1]
                    )
                    break
    
    return genes_candidates