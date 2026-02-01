#!/usr/bin/env python3
"""
Script to download sequences from NCBI GenBank
Reads sequence information from an input file
"""

import urllib.request
import urllib.parse
import argparse
import re
import time
import sys
import os

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
USER_AGENT = "Mozilla/5.0 (compatible; sequence-downloader/1.0)"

def fetch_genbank_record(accession):
    """
    Fetch the full GenBank record for an accession (needed to parse protein_id).
    
    Args:
        accession: GenBank accession number
        
    Returns:
        GenBank record as a string, or None on failure
    """
    params = {
        'db': 'nuccore',
        'id': accession,
        'rettype': 'gb',
        'retmode': 'text'
    }
    url = f"{EFETCH_URL}?{urllib.parse.urlencode(params)}"
    
    req = urllib.request.Request(url)
    req.add_header('User-Agent', USER_AGENT)
    with urllib.request.urlopen(req) as response:
        return response.read().decode('utf-8')

def parse_cds_location(genbank_record):
    """
    Parse CDS coordinates from the feature table of a GenBank record.
    Handles both simple (start..end) and multi-exon join(start1..end1,...) formats.
    The location may span multiple lines before the first qualifier line (starting with /).
    
    Args:
        genbank_record: Full GenBank record as a string
        
    Returns:
        List of (start, end) tuples using 1-based inclusive coordinates (as in GenBank),
        or None if no CDS feature is found
    """
    lines = genbank_record.split('\n')
    location_str = None
    in_cds = False

    for line in lines:
        if line.startswith('     CDS '):
            # First line of CDS feature; location begins at column 22
            in_cds = True
            location_str = line[21:].strip()
            # If there is no open join( or it is already closed, we are done
            if 'join(' not in location_str or ')' in location_str:
                break
        elif in_cds:
            # Continuation lines are indented; a qualifier line starts with /
            stripped = line.strip()
            if stripped.startswith('/'):
                break
            location_str += stripped
            if ')' in location_str:
                break

    if not location_str:
        return None

    # Strip join(...) wrapper if present
    location_str = location_str.replace('join(', '').rstrip(')')

    # Parse each start..end range
    result = []
    for part in location_str.split(','):
        part = part.strip()
        if '..' in part:
            start, end = part.split('..')
            result.append((int(start), int(end)))

    return result if result else None

def parse_origin_sequence(genbank_record):
    """
    Extract the full nucleotide sequence from the ORIGIN section of a GenBank record.
    The ORIGIN section contains numbered lines of sequence; this strips the numbers
    and whitespace to return a clean uppercase string.
    
    Args:
        genbank_record: Full GenBank record as a string
        
    Returns:
        Full nucleotide sequence as an uppercase string, or None if ORIGIN not found
    """
    in_origin = False
    parts = []
    for line in genbank_record.split('\n'):
        if line.startswith('ORIGIN'):
            in_origin = True
            continue
        if in_origin:
            if line.startswith('//'):
                break
            # Keep only alphabetic characters (strip position numbers and spaces)
            parts.append(re.sub(r'[^a-zA-Z]', '', line))

    return ''.join(parts).upper() if parts else None

def extract_cds_sequence(full_sequence, cds_ranges):
    """
    Extract the CDS nucleotides from the full sequence using GenBank coordinate ranges.
    For multi-exon genes the exon subsequences are concatenated in order.
    GenBank coordinates are 1-based and inclusive, so start is shifted by -1 for
    Python 0-based slicing while end is used as-is (making the slice inclusive).
    
    Args:
        full_sequence: The complete nucleotide sequence from the ORIGIN section
        cds_ranges: List of (start, end) tuples in 1-based inclusive coordinates
        
    Returns:
        The concatenated CDS nucleotide sequence
    """
    cds = ''
    for start, end in cds_ranges:
        cds += full_sequence[start - 1:end]
    return cds

def parse_protein_id(genbank_record):
    """
    Extract the protein_id from the CDS feature of a GenBank record.
    
    Args:
        genbank_record: Full GenBank record as a string
        
    Returns:
        Protein accession number (e.g., 'NP_002453.1'), or None if not found
    """
    match = re.search(r'/protein_id="([^"]+)"', genbank_record)
    if match:
        return match.group(1)
    return None

def fetch_protein_fasta(accession):
    """
    Fetch a protein FASTA sequence from NCBI.
    
    Args:
        accession: Protein accession number
        
    Returns:
        FASTA sequence as a string
    """
    params = {
        'db': 'protein',
        'id': accession,
        'rettype': 'fasta',
        'retmode': 'text'
    }
    url = f"{EFETCH_URL}?{urllib.parse.urlencode(params)}"

    req = urllib.request.Request(url)
    req.add_header('User-Agent', USER_AGENT)
    with urllib.request.urlopen(req) as response:
        return response.read().decode('utf-8')

def format_fasta(header, sequence):
    """
    Format a sequence as a FASTA string with 70-character lines.
    70 characters per line matches NCBI's own manual-download format.
    
    Args:
        header: FASTA header (without the leading >)
        sequence: Nucleotide or amino acid sequence string
        
    Returns:
        Complete FASTA-formatted string
    """
    fasta = f">{header}\n"
    for i in range(0, len(sequence), 70):
        fasta += sequence[i:i + 70] + "\n"
    return fasta

def parse_raw_fasta(raw_fasta):
    """
    Parse a raw FASTA string (as returned by NCBI) into its header and sequence.
    This is needed for the protein file: NCBI returns a pre-formatted FASTA that
    we must decompose so we can replace the header and re-wrap the sequence.
    
    Args:
        raw_fasta: Raw FASTA text (header line starting with >, followed by
                   sequence lines)
        
    Returns:
        Tuple of (header, sequence) where header is the text after > on the
        first line and sequence is the full amino-acid string with all
        line-breaks removed.
    """
    lines = raw_fasta.strip().split('\n')
    header = lines[0].lstrip('>')
    sequence = ''.join(lines[1:])
    return header, sequence

def protein_output_filename(output_filename):
    """
    Derive protein output filename from the nucleotide output filename.
    E.g., 'Human.fasta' -> 'Human_protein.fasta'
    
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
    E.g., 'human_mx1.fasta' -> 'human_mx1_stopcod.fasta'
    
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

def download_ncbi_sequence(accession, name, fasta_header, output_filename, output_dir, download_protein, rm_stop):
    """
    Download a coding sequence (CDS) from NCBI, and optionally its protein sequence.
    
    The script fetches the full GenBank record, then:
      1. Parses the CDS feature to get the exon coordinates
      2. Parses the ORIGIN section to get the full nucleotide sequence
      3. Extracts only the CDS portion (concatenating exons for multi-exon genes)
      4. Writes it as a FASTA file using fasta_header as the header line
    
    When rm_stop is True the last three nucleotides are checked for a stop codon.
    If one is found, two files are written: the original CDS (with stop codon)
    is saved with '_stopcod' in its filename, and the trimmed CDS (stop codon
    removed) is saved with the original filename. Both use fasta_header as header.
    
    When download_protein is True it also parses the protein_id from the same
    GenBank record and fetches the protein FASTA in a second request. The protein
    file uses fasta_header as its header.
    
    Args:
        accession: GenBank accession number
        name: Descriptive name for display
        fasta_header: Text to use after > in all FASTA headers for this sequence
        output_filename: Name for the nucleotide output file
        output_dir: Directory to save files
        download_protein: Whether to also download the protein sequence
        rm_stop: Whether to check for and remove stop codons
    """
    try:
        print(f"Downloading {name} ({accession})...")

        # Single fetch: the GenBank record contains everything we need
        print(f"  Fetching GenBank record...")
        genbank_record = fetch_genbank_record(accession)

        # --- Parse CDS coordinates ---
        cds_ranges = parse_cds_location(genbank_record)
        if not cds_ranges:
            print(f"  ✗ No CDS feature found in GenBank record for {accession}", file=sys.stderr)
            return False
        print(f"  CDS coordinates: {', '.join(f'{s}..{e}' for s, e in cds_ranges)}")

        # --- Parse the full sequence from the ORIGIN section ---
        full_sequence = parse_origin_sequence(genbank_record)
        if not full_sequence:
            print(f"  ✗ No ORIGIN sequence found in GenBank record for {accession}", file=sys.stderr)
            return False

        # --- Extract only the CDS nucleotides ---
        cds_sequence = extract_cds_sequence(full_sequence, cds_ranges)
        print(f"  Extracted CDS: {len(cds_sequence)} bp")

        # --- Write nucleotide CDS FASTA (with optional stop-codon handling) ---
        if rm_stop:
            trimmed_sequence, had_stop = check_and_remove_stop_codon(cds_sequence)
            if had_stop:
                stop_codon = cds_sequence[-3:]
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
            print(f"  ✓ Nucleotide CDS saved to {nuc_path}")

        # --- Protein: parse protein_id and fetch from NCBI protein db ---
        if download_protein:
            prot_id = parse_protein_id(genbank_record)
            if not prot_id:
                print(f"  ✗ No protein_id found in CDS feature for {accession}", file=sys.stderr)
                return False
            print(f"  Found protein_id: {prot_id}")

            time.sleep(0.5)
            prot_fasta_raw = fetch_protein_fasta(prot_id)
            _, prot_sequence = parse_raw_fasta(prot_fasta_raw)
            prot_path = os.path.join(output_dir, protein_output_filename(output_filename))
            with open(prot_path, 'w') as f:
                f.write(format_fasta(f"{fasta_header}", prot_sequence))
            print(f"  ✓ Protein saved to {prot_path}")

        return True

    except Exception as e:
        print(f"  ✗ Error downloading {name}: {str(e)}", file=sys.stderr)
        return False

def read_input_file(filename):
    """
    Read sequences from tab-delimited input file.
    
    Expected format (tab-delimited):
    accession_number    name    fasta_header    output_filename
    
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
                
                accession, name, fasta_header, output_filename = parts
                sequences.append({
                    'accession': accession.strip(),
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
        description='Download sequences from NCBI GenBank',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Input file format (tab-delimited, lines starting with # are comments):\n"
            "  accession_number\tname\tfasta_header\toutput_filename\n\n"
            "Examples:\n"
            "  python3 download_ncbi_sequences.py ncbi_sequences.tsv\n"
            "  python3 download_ncbi_sequences.py ncbi_sequences.tsv my_output_dir\n"
            "  python3 download_ncbi_sequences.py ncbi_sequences.tsv my_output_dir --protein\n"
            "  python3 download_ncbi_sequences.py ncbi_sequences.tsv my_output_dir --rmstop\n"
            "  python3 download_ncbi_sequences.py ncbi_sequences.tsv my_output_dir --protein --rmstop\n"
        )
    )
    parser.add_argument('input_file', help='Tab-delimited input file with sequence information')
    parser.add_argument('output_dir', nargs='?', default='.',
                        help='Output directory (created if it does not exist). Default: current directory')
    parser.add_argument('--protein', action='store_true',
                        help='Also download protein sequences. The script will parse the protein_id '
                             'from each GenBank record and download the corresponding protein FASTA. '
                             'Protein files are saved with "_protein" appended to the filename '
                             '(e.g., Human.fasta -> Human_protein.fasta)')
    parser.add_argument('--rmstop', action='store_true',
                        help='Check each nucleotide CDS for a stop codon (TAA, TAG, TGA) at the end '
                             'and remove it. When a stop codon is found two files are written: the '
                             'original sequence (with stop codon) is saved with "_stopcod" in its '
                             'filename, and the trimmed sequence (stop codon removed) is saved with '
                             'the original filename. Has no effect on protein sequences. '
                             '(e.g., Human_stopcod.fasta and Human.fasta)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    if args.output_dir != ".":
        try:
            os.makedirs(args.output_dir, exist_ok=True)
        except Exception as e:
            print(f"Error creating output directory '{args.output_dir}': {str(e)}", file=sys.stderr)
            sys.exit(1)
    
    print("=" * 60)
    print("NCBI Sequence Downloader")
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
        success = download_ncbi_sequence(
            seq['accession'],
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
        
        # Be nice to NCBI servers - wait between requests
        time.sleep(0.5)
    
    print()
    print("=" * 60)
    print(f"Download complete: {successful} successful, {failed} failed")
    print("=" * 60)

if __name__ == "__main__":
    main()
