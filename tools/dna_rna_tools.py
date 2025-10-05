def is_nucleic_acid(seq : str) -> bool:
    """
    Determines whether the sequence is nucleic acid. Case-sensitive.

    seq: str

    Returns bool.
    """
    valid_bases = "ATUCGatucg"
    has_t = any(base in "Tt" for base in seq)
    has_u = any(base in "Uu" for base in seq)
    if not all(base in valid_bases for base in seq):
        return False
    if has_t and has_u:
        return False
    else:
        return True


def transcribe(seq : str) -> str:
    """
    Transcribes the sequence by replacing "T" with "U" or vice versa.
    Case-sensitive.

    seq: str

    Returns str.
    """
    return seq.replace("T", "U").replace("t", "u")


def reverse(seq : str) -> str:
    """
    Reverses the sequence.

    seq: str
    Returns str.
    """
    return seq[::-1]


def complement(seq : str) -> str:
    """
    Provides a complementary sequence. Case-sensitive.

    seq: str

    Returns str.
    """
    dna_seq = any(base in "Tt" for base in seq)
    if dna_seq:
        complement_dict = str.maketrans("ATGCatgc", "TACGtacg")
    else:
        complement_dict = str.maketrans("AUGCaugc", "UACGuacg")
    return seq.translate(complement_dict)


def reverse_complement(seq : str) -> str:
    """
    Provides a reverse complementary sequence. Case-sensitive.

    seq: str

    Returns str.
    """
    return reverse(complement(seq))

