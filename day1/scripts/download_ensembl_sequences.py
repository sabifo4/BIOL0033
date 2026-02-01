#!/usr/bin/env python3
"""
Script to download coding sequences from Ensembl
Reads sequence information from an input file
Uses Ensembl REST API to retrieve coding sequences
"""

import urllib.request
import argparse
import json
import time
import sys
import os

USER_AGENT = "Mozilla/5.0 (compatible; sequence-downloader/1.0)"

def fetch_ensembl_sequence(transcript_id, seq_type):
    """
    Fetch a sequence from the Ensembl REST API.
    
    Args:
        transcript_id: Ensembl transcript ID (without version number)
        seq_type: 'cds' for nucleotide coding sequence, 'protein' for amino acid sequence
        
    Returns:
        Dictionary with 'id', 'desc', and 'seq' keys
    """
    url = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type={seq_type}"
    req = urllib.request.Request(url)
    req.add_header('Content-Type', 'application/json')
    req.add_header('User-Agent', USER_AGENT)
    
    with urllib.request.urlopen(req) as response:
        return json.loads(response.read().decode('utf-8'))

def format_fasta(header, sequence):
    """
    Format a sequence as a FASTA string with 60-character lines.
    
    Args:
        header: FASTA header text (without the leading >)
        sequence: Nucleotide or amino acid sequence string
        
    Returns:
        FASTA-formatted string
    """
    fasta = f">{header}\n"
    for i in range(0, len(sequence), 60):
        fasta += sequence[i:i+60] + "\n"
    return fasta

def protein_output_filename(output_filename):
    """
    Derive protein output filename from the nucleotide output filename.
    E.g., 'macaque_mx1.fasta' -> 'macaque_mx1_protein.fasta'
    
    Args:
        output_filename: Original nucleotide output filename
        
    Returns:
        Protein output filename
    """
    base, ext = os.path.splitext(output_filename)
    return f"{base}_protein{ext}"

def stopcod_output_filename(output_filename):
    """
    Derive the 'with stop codon' filename by inserting '_stopcod' before the extension.
    E.g., 'macaque_mx1.fasta' -> 'macaque_mx1_stopcod.fasta'
    
    Args:
        output_filename: Original output filename
        
    Returns:
        Filename with '_stopcod' inserted before the extension
    """
    base, ext = os.path.splitext(output_filename)
    return f"{base}_stopcod{ext}"

def check_and_remove_stop_codon(sequence):
    """
    Check whether the last three nucleotides of a CDS are a stop codon
    (TAA, TAG, or TGA) and, if so, return the sequence with it removed.
    
    Args:
        sequence: Uppercase nucleotide sequence string
        
    Returns:
        Tuple of (trimmed_sequence, had_stop_codon).
        trimmed_sequence is the sequence with the stop codon removed if one was
        present, or the original sequence unchanged if none was found.
        had_stop_codon is True when a stop codon was detected and removed.
    """
    stop_codons = {'TAA', 'TAG', 'TGA'}
    last_three = sequence[-3:].upper()
    if last_three in stop_codons:
        return sequence[:-3], True
    return sequence, False

def download_ensembl_sequence(gene_id, name, fasta_header, output_filename, output_dir, download_protein, rm_stop):
    """
    Download coding sequence for a gene from Ensembl, and optionally its protein sequence.
    
    The script first looks up the gene to find its canonical transcript, then fetches
    the CDS sequence. When download_protein is True, it also fetches the protein
    sequence using the same transcript ID with type=protein.

    All output FASTA files use fasta_header as their header line.

    When rm_stop is True the last three nucleotides are checked for a stop codon.
    If one is found, two files are written: the original CDS (with stop codon)
    is saved with '_stopcod' in its filename, and the trimmed CDS (stop codon
    removed) is saved with the original filename.
    
    Args:
        gene_id: Ensembl gene ID (e.g., ENSMMUG00000015329)
        name: Descriptive name for display
        fasta_header: Text to use after > in all FASTA headers for this sequence
        output_filename: Name for the nucleotide output file
        output_dir: Directory to save files
        download_protein: Whether to also download the protein sequence
        rm_stop: Whether to check for and remove stop codons
    """
    try:
        print(f"Downloading {name} ({gene_id})...")
        
        # Look up the gene to find the canonical transcript
        lookup_url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1"
        req = urllib.request.Request(lookup_url)
        req.add_header('Content-Type', 'application/json')
        req.add_header('User-Agent', USER_AGENT)
        
        with urllib.request.urlopen(req) as response:
            gene_data = json.loads(response.read().decode('utf-8'))
        
        canonical_transcript = gene_data.get('canonical_transcript')
        if not canonical_transcript:
            print(f"  ✗ No canonical transcript found for {gene_id}", file=sys.stderr)
            return False
        
        # Strip version number (e.g., ENSMMUT00000042651.3 -> ENSMMUT00000042651)
        transcript_id_base = canonical_transcript.split('.')[0]
        print(f"  Found canonical transcript: {canonical_transcript} (using {transcript_id_base})")
        time.sleep(0.3)

        # --- Download CDS (nucleotide) sequence (with optional stop-codon handling) ---
        cds_data = fetch_ensembl_sequence(transcript_id_base, 'cds')
        cds_sequence = cds_data['seq']

        if rm_stop:
            trimmed_sequence, had_stop = check_and_remove_stop_codon(cds_sequence)
            if had_stop:
                stop_codon = cds_sequence[-3:].upper()
                print(f"  Stop codon found: {stop_codon}")

                # Original (with stop codon) -> _stopcod filename
                stopcod_path = os.path.join(output_dir, stopcod_output_filename(output_filename))
                with open(stopcod_path, 'w') as f:
                    f.write(format_fasta(fasta_header, cds_sequence))
                print(f"  ✓ Nucleotide CDS with stop codon saved to {stopcod_path}")

                # Trimmed (stop codon removed) -> original filename
                nuc_path = os.path.join(output_dir, output_filename)
                with open(nuc_path, 'w') as f:
                    f.write(format_fasta(fasta_header, trimmed_sequence))
                print(f"  ✓ Nucleotide CDS without stop codon saved to {nuc_path}")
            else:
                print(f"  No stop codon found at end of CDS")
                nuc_path = os.path.join(output_dir, output_filename)
                with open(nuc_path, 'w') as f:
                    f.write(format_fasta(fasta_header, cds_sequence))
                print(f"  ✓ Nucleotide CDS saved to {nuc_path}")
        else:
            nuc_path = os.path.join(output_dir, output_filename)
            with open(nuc_path, 'w') as f:
                f.write(format_fasta(fasta_header, cds_sequence))
            print(f"  ✓ Nucleotide saved to {nuc_path}")

        # Download protein sequence if requested
        if download_protein:
            time.sleep(0.3)
            prot_data = fetch_ensembl_sequence(transcript_id_base, 'protein')
            prot_path = os.path.join(output_dir, protein_output_filename(output_filename))
            with open(prot_path, 'w') as f:
                f.write(format_fasta(f"{fasta_header}", prot_data['seq']))
            print(f"  ✓ Protein saved to {prot_path}")

        return True
        
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"  ✗ Gene ID {gene_id} not found in Ensembl", file=sys.stderr)
        else:
            print(f"  ✗ HTTP Error {e.code}: {str(e)}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"  ✗ Error downloading {name}: {str(e)}", file=sys.stderr)
        return False

def read_input_file(filename):
    """
    Read sequences from tab-delimited input file.
    
    Expected format (tab-delimited):
    gene_id    name    fasta_header    output_filename
    
    Args:
        filename: Path to input file
        
    Returns:
        List of dictionaries with sequence information
    """
    sequences = []
    
    try:
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                # Skip empty lines and comments
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Split by tab
                parts = line.split('\t')
                
                if len(parts) != 4:
                    print(f"Warning: Line {line_num} has {len(parts)} columns (expected 4), skipping", 
                          file=sys.stderr)
                    continue
                
                gene_id, name, fasta_header, output_filename = parts
                sequences.append({
                    'gene_id': gene_id.strip(),
                    'name': name.strip(),
                    'fasta_header': fasta_header.strip(),
                    'output_filename': output_filename.strip()
                })
        
        return sequences
        
    except FileNotFoundError:
        print(f"Error: Input file '{filename}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    """Main function to download all sequences"""
    parser = argparse.ArgumentParser(
        description='Download coding sequences from Ensembl',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Input file format (tab-delimited, lines starting with # are comments):\n"
            "  gene_id\tname\tfasta_header\toutput_filename\n\n"
            "Examples:\n"
            "  python3 download_ensembl_sequences.py ensembl_sequences.tsv\n"
            "  python3 download_ensembl_sequences.py ensembl_sequences.tsv my_output_dir\n"
            "  python3 download_ensembl_sequences.py ensembl_sequences.tsv my_output_dir --protein\n"
            "  python3 download_ensembl_sequences.py ensembl_sequences.tsv my_output_dir --rmstop\n"
            "  python3 download_ensembl_sequences.py ensembl_sequences.tsv my_output_dir --protein --rmstop\n"
        )
    )
    parser.add_argument('input_file', help='Tab-delimited input file with sequence information')
    parser.add_argument('output_dir', nargs='?', default='.',
                        help='Output directory (created if it does not exist). Default: current directory')
    parser.add_argument('--protein', action='store_true',
                        help='Also download protein sequences using the same canonical transcript. '
                             'Protein files are saved with "_protein" appended to the filename '
                             '(e.g., macaque_mx1.fasta -> macaque_mx1_protein.fasta)')
    parser.add_argument('--rmstop', action='store_true',
                        help='Check each nucleotide CDS for a stop codon (TAA, TAG, TGA) at the end '
                             'and remove it. When a stop codon is found two files are written: the '
                             'original sequence (with stop codon) is saved with "_stopcod" in its '
                             'filename, and the trimmed sequence (stop codon removed) is saved with '
                             'the original filename. Has no effect on protein sequences. '
                             '(e.g., macaque_mx1_stopcod.fasta and macaque_mx1.fasta)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    if args.output_dir != ".":
        try:
            os.makedirs(args.output_dir, exist_ok=True)
        except Exception as e:
            print(f"Error creating output directory '{args.output_dir}': {str(e)}", file=sys.stderr)
            sys.exit(1)
    
    print("=" * 60)
    print("Ensembl Sequence Downloader")
    print("=" * 60)
    print(f"Reading sequences from: {args.input_file}")
    print(f"Output directory:       {args.output_dir}")
    print(f"Download protein:       {'Yes' if args.protein else 'No'}")
    print(f"Remove stop codons:     {'Yes' if args.rmstop else 'No'}")
    print()
    
    # Read sequences from input file
    sequences = read_input_file(args.input_file)
    
    if not sequences:
        print("No sequences found in input file", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(sequences)} sequence(s) to download\n")
    
    successful = 0
    failed = 0
    
    for seq in sequences:
        success = download_ensembl_sequence(
            seq['gene_id'],
            seq['name'],
            seq['fasta_header'],
            seq['output_filename'],
            args.output_dir,
            args.protein,
            args.rmstop
        )
        
        if success:
            successful += 1
        else:
            failed += 1
        
        # Be nice to Ensembl servers - wait between requests
        time.sleep(0.5)
    
    print()
    print("=" * 60)
    print(f"Download complete: {successful} successful, {failed} failed")
    print("=" * 60)

if __name__ == "__main__":
    main()
