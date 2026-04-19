from abc import ABC, abstractmethod
from typing import Union, Tuple, Dict, Any
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord

class BiologicalSequence(ABC):
    """
    Abstract class for biological sequences.
    """
    def __init__(self, sequence: str):
        self._sequence = sequence
        if not self.check_alphabet():
            raise ValueError(f"Invalid characters in {self.__class__.__name__} sequence")
    
    def __len__(self) -> int:
        return len(self._sequence)
    
    def __getitem__(self, index: Union[int, slice]) -> Union[str, 'BiologicalSequence']:
        """
        Getting an element by index or slice.
        """
        result = self._sequence[index]
        if isinstance(index, slice):
            return self.__class__(result)
        return result
    
    def __str__(self) -> str:
        return self._sequence
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__} ({self._sequence})"
    
    @abstractmethod
    def check_alphabet(self) -> bool:
        """
        Checking whether the used alphabet is correct
        """
        pass

class NucleicAcidSequence(BiologicalSequence):
    """
    A class for nucleic acids (DNA and RNA)
    """
    _DNA_COMP_TABLE = str.maketrans("ATGCatgc", "TACGtacg")
    _RNA_COMP_TABLE = str.maketrans("AUGCaugc", "UACGuacg")

    @property
    @abstractmethod
    def _complement_table(self):
        pass

    @property
    @abstractmethod
    def _alphabet(self) -> str:
        pass

    def check_alphabet(self) -> bool:
        """
        Checking whether all symbols are in the correct alphabet
        for nucleic acids.
        """
        return all(base in self._alphabet for base in self._sequence)
    
    def complement(self) -> 'NucleicAcidSequence':
        """
        Returns a complementary sequence
        """
        complemented = self._sequence.translate(self._complement_table)
        return self.__class__(complemented)
    
    def reverse_complement(self) -> 'NucleicAcidSequence':
        """
        Returns a reverse complementary sequence
        """
        return self.reverse().complement()
    
    def reverse(self) -> 'NucleicAcidSequence':
        """
        Returns a reversed sequence
        """
        return self.__class__(self._sequence[::-1])

class DNASequence(NucleicAcidSequence):
    """
    A class for DNA sequences
    """
    _alphabet = "ATGCatgc"
    _complement_table = NucleicAcidSequence._DNA_COMP_TABLE

    def transcribe(self) -> 'RNASequence':
        """
        Transcribes DNA into RNA
        """
        transcribed = self._sequence.replace("T", "U").replace("t", "u")
        return RNASequence(transcribed)

class RNASequence(NucleicAcidSequence):
    """
    A class for RNA sequences
    """
    _alphabet = "AUGCaugc"
    _complement_table = NucleicAcidSequence._RNA_COMP_TABLE

    def reverse_transcribe(self) -> 'DNASequence':
        """
        Reverse transcribes RNA into DNA
        """
        reverse_transcribed = self._sequence.replace("U", "T").replace("u", "t")
        return DNASequence(reverse_transcribed)

class AminoAcidSequence(BiologicalSequence):
    """
    A class for amino acid sequences
    """
    _alphabet = "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy"

    def check_alphabet(self) -> bool:
        """
        Checking whether all symbols are in the correct alphabet
        for amino acids.
        """
        return all(aa in self._alphabet for aa in self._sequence)
    
    def calculate_molecular_weight(self) -> float:
        """
        Calculates approximate protein molecular mass
        """
        weights = {
            'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
            'Q': 146.15, 'E': 147.13, 'G': 75.07, 'H': 155.16, 'I': 131.17,
            'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
            'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
        }
        total_weight = 0.0
        for aa in self._sequence.upper():
            if aa in weights:
                total_weight += weights[aa]
            else:
                total_weight += 110.0 # Average value for non-canonic AAs
        water_weight = 18.015 * (len(self._sequence) - 1)
        return total_weight - water_weight



def define_bounds(bounds: Union[Tuple[float, float], int, float]) -> Tuple[float, float]:
    """
    Defines bounds as a min and a max depending on the type

    Args:
    bounds: int / float (interpreted as the upper bound, the lower bound is 0) or tuple (min, max)

    Returns tuple.
    """
    if isinstance(bounds, (int, float)):
        return (0.0, float(bounds))
    elif isinstance(bounds, tuple) and len(bounds) == 2:
        return (float(bounds[0]), float(bounds[1]))
    return bounds


def calculate_quality(record: SeqRecord) -> float:
    """
    Calculates average sequence quality.

    Args:
    record: SeqRecord object with phred_quality annotation

    Returns avg quality
    """
    if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
        qualities = record.letter_annotations['phred_quality']
        return sum(qualities) / len(qualities)
    return 0.0


def filter_fastq(input_file: str, output_file: str,
                 gc_bounds: Union[Tuple[float, float], float, int] = (0, 100),
                 length_bounds: Union[Tuple[int, int], float, int] = (0, 2**32),
                 quality_threshold: Union[float, int] = 0,
                 log_file: str = 'filter_fastq.log') -> Dict[str, int]:
    """
    Filters FastQ sequences based on GC content, length, and quality, 
    with user-defined bounds and thresholds.

    Args:
    input_file: path to input fastq file
    output_file: path to output fastq file
    gc_bounds: int / float / tuple (default (0, 100))
    length_bounds: int / float / tuple (default (0, 2**32))
    quality_threshold: int / float (default 0)

    Returns dict.
    """
    gc_min, gc_max = define_bounds(gc_bounds)
    len_min, len_max = define_bounds(length_bounds)
    stats = {
        'total': 0,
        'passed': 0,
        'filtered_gc': 0,
        'filtered_length': 0,
        'filtered_quality': 0
    }
    with open(input_file, 'r') as in_handle:
        with open(output_file, 'w') as out_handle:
            for record in SeqIO.parse(in_handle, 'fastq'):
                stats['total'] += 1
                seq_len = len(record.seq)
                if not (len_min <= seq_len <= len_max):
                    stats['filtered_length'] += 1
                    continue
                gc_content = gc_fraction(record.seq) * 100
                if not (gc_min <= gc_content <= gc_max):
                    stats['filtered_gc'] += 1
                    continue
                avg_quality = calculate_quality(record)
                if avg_quality < quality_threshold:
                    stats['filtered_quality'] += 1
                    continue
                SeqIO.write(record, out_handle, 'fastq')
                stats['passed'] += 1
    return stats