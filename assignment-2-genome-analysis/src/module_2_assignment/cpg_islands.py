def cpg_islands(  sequence: str, window_size: int = 200, 
                gc_threshold: float = 0.5, 
                obs_exp_threshold: float = 0.6):

    sequence = sequence.upper()
    seq_length = len(sequence)
    results = []

    if seq_length < window_size:
        return results

    # Calculate counts for the first window
    window_seq = sequence[:window_size]

    c_count = window_seq.count("C")
    g_count = window_seq.count("G")
    cg_count = window_seq.count("CG")

    # Slide the window across the sequence
    for i in range(seq_length - window_size + 1):

        # Calculate CpG statistics
        gc_count = c_count + g_count
        gc_content_value = gc_count / window_size

        expected_cg = (c_count * g_count) / window_size

        obs_exp_ratio = (
            cg_count / expected_cg
            if expected_cg > 0
            else 0
        )

        if ( gc_content_value > gc_threshold and obs_exp_ratio > obs_exp_threshold):
            
            results.append( ( i + 1,  
                             i + window_size, 
                             obs_exp_ratio, 
                             gc_content_value ) )

        # Update counts for next window
        if i < seq_length - window_size:

            outgoing_base = sequence[i]
            incoming_base = sequence[i + window_size]

            # Remove outgoing base
            if outgoing_base == "C":
                c_count -= 1
            elif outgoing_base == "G":
                g_count -= 1

            # Add incoming base
            if incoming_base == "C":
                c_count += 1
            elif incoming_base == "G":
                g_count += 1

            # Remove CG pair leaving the window
            if sequence[i:i + 2] == "CG":
                cg_count -= 1

            # Add CG pair entering the window
            new_pair_start = i + window_size - 1

            if sequence[new_pair_start:new_pair_start + 2] == "CG":
                cg_count += 1

    return results