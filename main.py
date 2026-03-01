from tools.dna_rna_tools import (is_nucleic_acid, transcribe, reverse, 
                           complement, reverse_complement)
from tools.fastqc_tools import (define, check_gc, check_length, check_quality)


def run_dna_rna_tools(*args : list[str]) -> list[str] | str:
    """
    Accepts nucleic acid sequences and performs procedures with them,
    such as determining whether it's a nucleic acid sequence,
    returning a transcribed, reverse, complementary, or reverse complementary sequences.
    Case-sensitive.

    args: list[str]

    Returns list[str] / str.
    Raises exception if less then 2 args are provided, if an unknown procedure is used,
    or if one or more of the sequences are not nucleic acids
    (and the procedure is not "is_nucleic_acid").
    """
    if len(args) < 2:
        raise ValueError("At least 1 sequence and procedure are expected")
    sequences = args[:-1]
    procedure = args[-1]
    results = []
    for seq in sequences:
        if procedure == "is_nucleic_acid":
            results.append(is_nucleic_acid(seq))
        elif not is_nucleic_acid(seq):
            raise ValueError("Only nucleic acid sequences are accepted")
        else:
            if procedure == "transcribe":
                results.append(transcribe(seq))
            elif procedure == "reverse":
                results.append(reverse(seq))
            elif procedure == "complement":
                results.append(complement(seq))
            elif procedure == "reverse_complement":
                results.append(reverse_complement(seq))
            else:
                raise ValueError("Unknown procedure")
    return results[0] if len(results) == 1 else results


def filter_fastq(seqs : dict, gc_bounds : tuple | int | float = (0, 100), 
                 length_bounds : tuple | int | float = (0, 2**32),
                 quality_threshold : int | float = 0) -> dict:
    """
    Filters FastQ sequences based on GC content, length, and quality, 
    with user-defined bounds and thresholds.

    seqs: dict {name: str: (sequence: str, quality_str: str)}
    gc_bounds: int / float / tuple (default (0, 100))
    length_bounds: int / float / tuple (default (0, 2**32))
    quality_threshold: int / float (default 0)

    Returns dict.
    """
    gc_min, gc_max = define(gc_bounds)
    len_min, len_max = define(length_bounds)
    filtered_seqs = {}
    for name, (sequence, quality_str) in seqs.items():
        if (check_gc(sequence, gc_min, gc_max) and
           check_length(sequence, len_min, len_max) and
           check_quality(quality_str, quality_threshold)):
            filtered_seqs[name] = (sequence, quality_str)
    return filtered_seqs
