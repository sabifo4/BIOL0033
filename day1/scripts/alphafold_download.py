#!/usr/bin/env python3
"""
Script to download metadata and PDB files from alphafold for uniprot identifiers
Ouputs .json files named <uniprot_id>.json in the specified output directory.

"""

import os
import pandas as pd
import numpy as np
import argparse
import sys
import json
import requests



def download_alphafold_metadata(uniprot_id, species_name, model_metadata_dir):
    # AlphaFold DB API URL for metadata
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        
    # Output file path
    output_file = os.path.join(model_metadata_dir, f"metadata_{species_name}_{uniprot_id}.json")
    
    try:
        print(f"Downloading metadata for {uniprot_id} ({species_name})...")
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad status codes
        
        # Convert response to JSON format
        metadata = response.json()
        
        # Save the JSON file
        with open(output_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"  ✓ Saved metadata to {output_file}")
        
    except requests.exceptions.RequestException as e:
        print(f"  ✗ Failed to download metadata for {uniprot_id}: {e}")
    except json.JSONDecodeError as e:
        print(f"  ✗ Failed to decode JSON for {uniprot_id}: {e}")

    return

def download_pdb(uniprot_id, species_name, model_metadata_dir, pdb_dir):
    output_file = os.path.join(pdb_dir, f"{species_name}_{uniprot_id}.pdb")

    metadata = json.load(open(os.path.join(model_metadata_dir, f"metadata_{species_name}_{uniprot_id}.json")))
    pdb_url = metadata[0]['pdbUrl']

    try:
        print(f"Downloading pdb for {uniprot_id} ({species_name})...")
        response = requests.get(pdb_url)
        response.raise_for_status()  # Raise an error for bad status codes
        
        pdb_data = response.text
        
        # Save the PDB file
        with open(output_file, 'w') as f:
            f.write(pdb_data)
        print(f"  ✓ Saved PDB to {output_file}")
        
    except requests.exceptions.RequestException as e:
        print(f"  ✗ Failed to download PDB for {uniprot_id}: {e}")

def main():
    """Main function to download structures and structure metadata from AlphaFold DB."""
    parser = argparse.ArgumentParser(
        description='Download structures and structure metadata from AlphaFold DB',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Input file format (csv):\n"	
            "  species_common_name\tgene_name\tdatabase\taccession\tuniprot\n"
            "Examples:\n"
            "  python3 alphafold_download.py protein_metadata.csv\n"
            "  python3 alphafold_download.py protein_metadata.csv my_output_dir\n"
        )
    )
    parser.add_argument('input_file', help='Tab-delimited input file with sequence information')
    parser.add_argument('output_dir', nargs='?', default='.',
                        help='Output directory (created if it does not exist). Default: current directory')
    # parser.add_argument('--protein', action='store_true',
    #                     help='Also download protein sequences. The script will parse the protein_id '
    #                          'from each GenBank record and download the corresponding protein FASTA. '
    #                          'Protein files are saved with "_protein" appended to the filename '
    #                          '(e.g., Human.fasta -> Human_protein.fasta)')
    # parser.add_argument('--rmstop', action='store_true',
    #                     help='Check each nucleotide CDS for a stop codon (TAA, TAG, TGA) at the end '
    #                          'and remove it. When a stop codon is found two files are written: the '
    #                          'original sequence (with stop codon) is saved with "_stopcod" in its '
    #                          'filename, and the trimmed sequence (stop codon removed) is saved with '
    #                          'the original filename. Has no effect on protein sequences. '
    #                          '(e.g., Human_stopcod.fasta and Human.fasta)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    if args.output_dir != ".":
        try:
            os.makedirs(args.output_dir, exist_ok=True)
        except Exception as e:
            print(f"Error creating output directory '{args.output_dir}': {str(e)}", file=sys.stderr)
            sys.exit(1)
    
    print("=" * 60)
    print("AlphaFold PDB Downloader")
    print("=" * 60)
    print(f"Reading sequences from: {args.input_file}")
    print(f"Output directory:       {args.output_dir}")
    print()
    
    # Create model_metadata directory if it doesn't exist
    model_metadata_dir = os.path.join(args.output_dir, "model_metadata")
    try:
        os.makedirs(model_metadata_dir, exist_ok=True)
    except Exception as e:
        print(f"Error creating model_metadata directory '{model_metadata_dir}': {str(e)}", file=sys.stderr)
        sys.exit(1)

    # Create pdb directory
    pdb_dir = os.path.join(args.output_dir, "pdb")
    try:
        os.makedirs(pdb_dir, exist_ok=True)
    except Exception as e:
        print(f"Error creating pdb directory '{pdb_dir}': {str(e)}", file=sys.stderr)
        sys.exit(1)

    # successful = 0
    # failed = 0
    
    protein_metadata = pd.read_csv(args.input_file, dtype=str)
    
    for idx, row in protein_metadata.iterrows():
        uniprot_id = row['uniprot'].strip()  # Remove any whitespace
        species_name = row['species_common_name']
        
        download_alphafold_metadata(uniprot_id, species_name, model_metadata_dir)
        download_pdb(uniprot_id, species_name, model_metadata_dir, pdb_dir)
    
    print(f"\nMetadata and PDB download complete!")
    print(f"\nMetadata Files saved in: {model_metadata_dir}")
    print(f"\nPDB Files saved in: {pdb_dir}")

    # for seq in sequences:
    #     success = download_ncbi_sequence(
    #         seq['accession'],
    #         seq['name'],
    #         seq['fasta_header'],
    #         seq['output_filename'],
    #         args.output_dir,
    #         args.protein,
    #         args.rmstop
    #     )
        
    #     if success:
    #         successful += 1
    #     else:
    #         failed += 1
        
    #     # Be nice to NCBI servers - wait between requests
    #     time.sleep(0.5)
    
    # print()
    # print("=" * 60)
    # print(f"Download complete: {successful} successful, {failed} failed")
    # print("=" * 60)

if __name__ == "__main__":
    main()