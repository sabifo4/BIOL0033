#!/usr/bin/env python3

import argparse
import os
import subprocess
from Bio import SeqIO

"""
Script to download the translated pig sequence for M65087 from NCBI and replace the translated AB164037 sequence in the ../raw_data/data1/unaln_nuc.fasta file with it.
"""

def main():
    """Takes input unaln_nuc.fasta and outputs a new file with the translated pig sequence for M65087 replacing the translated AB164037 sequence with the new filename <input_file>_pigseqswap.fasta in the <output_directory>."""
    parser = argparse.ArgumentParser(
        description='download and swap pig sequence',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Input file:  fasta listing unaligned nucleotide sequences\n"	
            "Example:\n"
            "  python3 unaln_nuc.fasta unaln_nuc_pigseqswap\n"
        )
    )
    parser.add_argument('input_file', help='input unaligned nucleotide fasta file')
    parser.add_argument('output_directory', help='directory to save output file')

    args = parser.parse_args()

    print(args.output_directory)
    print(os.path.join(args.output_directory, "Pig.fasta"))
    
    #First make a small file as input to the download_ncbi_sequence.py script to download the pig sequence for M65087
    with open(os.path.join(args.output_directory, 'pig_fasta_metadata.txt'), 'w') as f:    
            f.write("# NCBI Mx1 sequences\n")
            f.write("# Format: accession_number	name	fasta_header	output_filename\n")
            f.write("M65087	Sus_scrofa	Pig_Mx	Pig.fasta\n")
    
    #Download the pig sequence for M65087 using the download_ncbi_sequence.py script
    subprocess.run(["python3", "../../scripts/download_ncbi_sequences.py", os.path.join(args.output_directory, "pig_fasta_metadata.txt"), args.output_directory, "--protein", "--rmstop"], cwd='.')
    #subprocess.run(["python3", "subprocess_test.py"], cwd='../../scripts/')
    #subprocess.call(["python3", "~/my_session/day1/scripts/download_fasta.py", os.path.join(args.output_directory, "pig_fasta_metadata.txt"), args.output_directory, "--protein", "--rmstop"])
    

    #Replace the pig sequence and write to new output file
    output_fasta = SeqIO.parse(os.path.join(args.input_file), "fasta")
    pig_fasta = SeqIO.parse(os.path.join(args.output_directory, "Pig.fasta"), "fasta")
    pig_record = next(pig_fasta)

    with open(os.path.join(args.output_directory, os.path.basename(args.input_file).split('.')[0] + "_pigseqswap.fasta"), 'w') as f_out:
        for record in output_fasta:
            f_out.write('>' + record.id + '\n')
            if record.id == "Pig_Mx":
                f_out.write(str(pig_record.seq) + '\n')
            else:
                f_out.write(str(record.seq) + '\n')
            


if __name__ == "__main__":
    main()