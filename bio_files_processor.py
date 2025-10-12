import os
import re

def convert_multiline_fasta_to_one_line(input_fasta : str, output_fasta : str | None = None) -> None:
    """
    Converts multiline sequences in FASTA files to a single line.

    input_fasta: str
    output_fasta: str / None

    Returns None.
    """
    if output_fasta is None:
        base_input_name = os.path.basename(input_fasta)
        name, my_extension = os.path.splitext(base_input_name)
        output_fasta = f"{name}_oneline{my_extension}"
    if os.path.abspath(input_fasta) == os.path.abspath(output_fasta):
        raise ValueError("The input file cannot be overwritten.")
    if os.path.exists(output_fasta):
        raise FileExistsError("File already exists.")
    with open(input_fasta, 'r') as input_file, open(output_fasta, 'a') as output_file:
        header = None
        sequence = []
        for line in input_file:
            line = line.strip()
            if line.startswith('>'):
                header = line
                sequence = []
                if header is not None:
                    output_file.write(f"{header}\n")
                    output_file.write(f"{''.join(sequence)}\n")
            else:
                sequence.append(line)
        if header is not None:
            output_file.write(f"{header}\n")
            output_file.write(f"{''.join(sequence)}\n")


def parse_blast_output(input_file : str, output_file : str) -> None:
    """
    Parses Blast output to find the description of the best match.

    input_file: str
    output_file: str

    Returns None.
    """
    if os.path.abspath(input_file) == os.path.abspath(output_file):
        raise ValueError("The input file cannot be overwritten.")
    if os.path.exists(output_file):
        raise FileExistsError("File already exists.")
    best_matches = []
    with open(input_file, 'r') as input:
        lines = input_file.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "Sequences producing significant alignments:":
            i += 3 # Skipping the blank line and table header
            if i < len(lines):
                match = lines[i].strip()
                split_match = re.split("... |  ", match)
                description = split_match[0]
                best_matches.append(description)
        else:
            i += 1
    best_matches.sort()
    with open(output_file, 'w') as output:
        for match in best_matches:
            output.write(f"{match}\n")

