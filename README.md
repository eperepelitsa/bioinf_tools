# bioinf_tools
Bioinformatic tools for work with DNA/RNA/protein sequences and FastQ files.

## Installation
To run bioinf_tools, you need to clone it from the repository:
```
git clone https://github.com/eperepelitsa/bioinftools.git
cd biotools
```

## Usage
### Work with Biological Sequences
The module provides an object-oriented interface for working with different biological sequences.
#### Creating Sequence Objects
```
from bioinf_tools import DNASequence, RNASequence, AminoAcidSequence

# Create DNA sequence
dna = DNASequence("ATGCatgc")
print(dna) # ATGCatgc

# Create RNA sequence
rna = RNASequence("AUGCaugc")

# Create protein sequence
protein = AminoAcidSequence("MQARSRD")
```

#### Sequence Operations:
All sequence types support:
* Getting length: `len(dna)`
* Indexing and slicing: `dna[2], dna[1:4]`
* String representation: `str(dna)`

#### DNA Sequence Methods
```
dna = DNASequence("ATGCatgc")

# Complementary sequence
dna.complement() # TACGtacg

# Reverse complement
dna.reverse_complement() # gcatGCAT

# Transcribe DNA to RNA
dna.transcribe() # AUGCaugc
```

#### RNA Sequence Methods
```
rna = RNASequence("AUGCaugc")

# Complementary sequence
rna.complement() # UACGuacg

# Reverse complement
rna.reverse_complement() # gcauGCAU

# Transcribe RNA to DNA
rna.reverse_transcribe() # ATGCatgc
```

#### Protein Sequence Methods
```
protein = AminoAcidSequence("MQARSRD")

# Calculate approximate molecular weight (in Da)
protein.calculate_molecular_weight() # 862.95
```

### FastQ Filtration
Function `filter_fastq` filters FastQ sequences depending on GC content, length, and quality.
```
from bioing_tools import filter_fastq

stats = filter_fastq(
    input_file="sample.fastq"
    output_file="filtered_fastq",
    gc_bounds=(40,60),             # GC content between 40% and 60%
    length_bounds=(50,500),        # Sequence length between 50 and 500 bp
    quality_threshold=20           # Average quality score >= 20
)
print(stats)
# {'total': 1000, 'passed': 750, 'filtered_gc': 150,
# 'filtered_length': 50, 'filtered_quality': 50}
```

#### Parameter Options
##### GC Bounds (`gc_bounds`)
* As tuple: `gc_bounds=(40,60)` - min and max GC content
* As number: `gc_bounds=60` - max GC content (min = 0)

##### Length Bounds (`length_bounds`)
* As tuple: `length_bounds=(50,500)` - min and max length
* As number: `length_bounds=300` - max length (min = 0)

##### Quality threshold (`quality_threshold`)
* Single number: `quality_threshold=25` - min avg quality

## Requirements
* Python 3.6+
* Biopython

## Contacts
- Author: Elizaveta Perepelitsa
- Email: ross8marquise@gmail.com
