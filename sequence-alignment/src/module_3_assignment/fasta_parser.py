import warnings

# Considering valid characters for each type sa per the table above
VALID_DNA = set("ACGTMSWBDRHVYKN-")
VALID_RNA = set("ACGUMSWBDRHVYKN-")
VALID_PROTEIN = set("ABCDEFGHIKLMNPQRSTUVWYZX*-")

# Checking whether a sequence contains valid biological characters
def check_valid(sequence):

    # if sequence is empty, its invalid
    if not sequence:
        return False

    sequence = sequence.upper()

    # a sequence is valid, when all of its characters are present
    #  in atleast one of the 3 types above
    if all(char in VALID_DNA for char in sequence):
        return True

    if all(char in VALID_RNA for char in sequence):
        return True

    if all(char in VALID_PROTEIN for char in sequence):
        return True

    return False

# Classifying the type of sequence 
def classify(sequence):

    if not sequence:
        return "Unknown"

    sequence = sequence.upper()

    if all(char in VALID_DNA for char in sequence):
        return "DNA"

    elif all(char in VALID_RNA for char in sequence):
        return "RNA"

    elif all(char in VALID_PROTEIN for char in sequence):
        return "Protein"

    else:
        return "Unknown"

def read_fasta(file_path):

    description = None
    sequence = ""

    with open(file_path, "r") as f:

        for line_number, line in enumerate(f, start=1):

            # removing whitespace from end of line
            line = line.rstrip()

            # Ignoring blank lines
            if not line:
                continue

            # Reading the header line
            if line.startswith(">"):

                # If we already have a record, we process it first before starting 
                # a new one
                if description is not None:

                    if not check_valid(sequence):
                        warnings.warn(
                            f"Invalid sequence found for "
                            f"'{description}'. Skipping it."
                        )
                    else:
                        seq_type = classify(sequence)
                        if seq_type == "Unknown":
                            warnings.warn(
                                f"Unknown sequence type found for "
                                f"'{description}'. Skipping it."
                            )
                        else:
                            yield (description, sequence)

                # Starting a new record , everything after ">" is description
                description = line[1:]

                # reseting sequence for new record
                sequence = ""

            # Reading the sequence lines
            else:

                # Sequence data before first header is invalid
                if description is None:
                    warnings.warn(
                        f"Line {line_number}: sequence data found "
                        f"before first FASTA header. Skipping line."
                    )
                    continue

                # Adding line to current sequence
                sequence += line.upper()

        # Processing the final fasta record
        if description is not None:

            if not check_valid(sequence):
                warnings.warn(
                    f"Invalid sequence found for "
                    f"'{description}'. Skipping it."
                )
            else:
                seq_type = classify(sequence)
                if seq_type == "Unknown":
                    warnings.warn(
                        f"Unknown sequence type found for "
                        f"'{description}'. Skipping it."
                    )
                else:
                    yield (description, sequence)
