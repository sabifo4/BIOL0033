#!/usr/bin/env python3

import argparse

"""
Script to remove header and footer from USalign
"""

def main():
    """Takes input file and outputs a parsed alignment file in basic fasta format (No header and footer)."""
    parser = argparse.ArgumentParser(
        description='Convert USalign fasta to remove header',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Input file:  fasta with header from USalign\n"	
            "Example:\n"
            "  python3 usalign_fasta_parse.py input.fasta output.fasta\n"
        )
    )
    parser.add_argument('input_file', help='input alignment from USalign')
    parser.add_argument('output_file', help='filename')

    args = parser.parse_args()

    with open(args.input_file, 'r') as f_in: 
        with open(args.output_file, 'w') as f_out:
            for line in f_in: 
                if line[0]==">":
                    f_out.write(line)
                    seq_line = next(f_in)
                    f_out.write(seq_line)
                
                

if __name__ == "__main__":
    main()