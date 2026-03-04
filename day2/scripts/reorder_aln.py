#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO

"""
Script to reorder the sequences in the input fasta to match the order of the sequences in the reference fasta fastaa_aln.fasta file. The output file will be saved as <input_file>_reordered.fasta in the <output_directory>.
"""

def main():
    parser = argparse.ArgumentParser(
        description='reorder sequences in input fasta to match order of sequences in reference fasta',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "input_fasta: fasta file with sequences out of order\n"	
            "reference_fasta: fasta file with sequences in correct order\n"
            "output_directory: directory to save output file\n"
            "Example:\n"
            "  python3 reorder_aln.py aa_seq/aa_aln_seq.fasta aa/aa_aln.fasta aa_seq\n"
        )
    )
    parser.add_argument('input_file', help='fasta file with sequences out of order')
    parser.add_argument('reference_file', help='reference fasta file with sequences in correct order')
    parser.add_argument('output_directory', help='directory to save output file')

    args = parser.parse_args()
    
    #Read in the input fasta as a dictionary
    input_fasta = SeqIO.to_dict(SeqIO.parse(os.path.join(args.input_file), "fasta"))

    #Read in the reference fasta to get the order of the sequences
    reference_fasta = SeqIO.parse(os.path.join(args.reference_file), "fasta")

    #get base filename without extension and directory for output file
    base_output_filename = os.path.basename(args.input_file).split('/')[-1].split('.')[0]
    output_filename = os.path.join(args.output_directory, base_output_filename + "_reordered.fasta")

    #Write the sequences in the order of the reference fasta to the output file
    with open(output_filename, 'w') as f_out:
        for record in reference_fasta:
            if record.id in input_fasta:
                f_out.write('>' + record.id + '\n')
                f_out.write(str(input_fasta[record.id].seq) + '\n')
            else:
                print(f"Warning: {record.id} not found in input fasta. Skipping this sequence.")
    
    print(f"Reordered fasta file saved to {output_filename}")

if __name__ == "__main__":
    main()