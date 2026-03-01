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
    return avg_quality >= quality_threshold

