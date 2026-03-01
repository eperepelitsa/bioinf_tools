# bioinf_tools
Bioinformatic tools for work with DNA/RNA sequences and FastQ files.

## Installation
To run bioinf_tools, you need to clone it from the repository:
```
git clone https://github.com/eperepelitsa/bioinftools.git
cd biotools
```

## Usage
### Work with DNA/RNA sequences
Function `run_dna_rna_tools` accepts 1 or more sequences and 1 of the following procedures:
- Verifying the nucleic acid sequence
  ```
  run_dna_rna_tools('TTUU', 'is_nucleic_acid') # False
  ```
- DNA -> RNA transcription
  ```
  run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
  ```
- Providing the complementary sequence
  ```
  run_dna_rna_tools('AtG', 'complement') # 'TaC'
  ```
- Reversing the sequence
  ```
  run_dna_rna_tools('ATG', 'reverse') # 'GTA'
  ```
- Providing the reverse complementary sequence
  ```
  run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
  ```

### Filtration of FastQ sequences
Function `filter_fastq` filters FastQ sequences depending on:
- GC content
  ```
  filter_fastq(sequences, gc_bounds = (40, 60))
  ```
- length
  ```
  filter_fastq(sequences, length_bounds = 100)
  ```
- quality
  ```
  filter_fastq(sequences, quality_threshold = 20)
  ```

  ## Contacts
  - Author: Elizaveta Perepelitsa
  - Email: ross8marquise@gmail.com
