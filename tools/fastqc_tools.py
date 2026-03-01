def define(bounds : tuple | int | float) -> bool:
    """
    Defines bounds as a min and a max depending on the type

    Arguments:
    bounds: int / float / tuple

    Returns tuple.
    """
    if type(bounds) in (int, float):
        return (0, bounds)
    elif type(bounds) == tuple and len(bounds) == 2:
        return bounds


def check_gc(sequence : str, gc_min : int | float, gc_max : int | float) -> bool:
    """
    Calculates GC content and checks whether it's acceptable depending on bounds.

    Arguments:
    sequence: str
    gc_min: int / float
    gc_max: int / float

    Returns bool.
    """
    gc_count = sequence.upper().count("G") + sequence.upper().count("C")
    gc_content = (gc_count/len(sequence)) * 100
    if len(sequence) == 0:
        gc_content = 0
    return gc_min <= gc_content <= gc_max


def check_length(sequence : str, len_min : int | float, len_max : int | float) -> bool:
    """
    Calculates sequence length and checks whether it's acceptable depending on bounds.

    Arguments:
    sequence: str
    len_min: int
    len_max: int

    Returns bool.
    """
    return len_min <= len(sequence) <= len_max


def check_quality(quality_str : str, quality_threshold : int | float) -> bool:
    """
    Calculates average sequence quality and checks
    whether it's acceptable depending on the threshold.

    Arguments:
    quality_str: str
    quality_threshold: int / float

    Returns bool.
    """
    quality_scores = [ord(char) - 33 for char in quality_str]
    avg_quality = sum(quality_scores)/len(quality_scores)
    if len(quality_scores) == 0:
        avg_quality = 0
    return avg_quality >= quality_threshold


def read_fastq(input_fastq : str) -> dict:
    """
    Reads FastQ files line by line to find
    the sequence ID, the sequence and the quality string.

    Arguments:
    input_fastq: str

    Returns dict.
    """
    seqs = {}
    with open(input_fastq, 'r') as file:
        while True:
            seq_id = file.readline().strip()
            if not seq_id:
                break  # End of the file
            sequence = file.readline().strip()
            plus_line = file.readline().strip()
            quality_str = file.readline().strip()
            if not seq_id.startswith('@') or not plus_line.startswith('+'):
                raise ValueError("Only FastQ files are accepted")
            seqs[seq_id] = (sequence, quality_str)
    return seqs


def write_fastq(output_fastq : str, filtered_seqs : dict) -> None:
    """
    Writes the filtered results in a separate file.

    Arguments:
    output_fastq: str
    filtered_seqs: dict

    Returns None.
    """
    with open(output_fastq, 'w') as file:
        for seq_id, sequence, plus_line, quality_str in filtered_seqs:
            file.write(f"{seq_id}\n")
            file.write(f"{sequence}\n")
            file.write(f"{plus_line}\n")
            file.write(f"{quality_str}\n")
