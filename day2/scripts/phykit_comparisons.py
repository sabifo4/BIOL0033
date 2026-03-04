#!/usr/bin/env python3

#Code was written by Dr. Ben Heineike, with assitance from Claude Haiku 4.5.

import os
import argparse
import subprocess
import csv
from Bio import SeqIO
from pathlib import Path

"""
Script to run metrics from phykit on alignments and trees from a list of input alignments, a selection of phykit metrics, and tree files. The output will be saved as output_filename.txt.
"""

def parse_input_files(input_files):
    """
    Parse the input_files file and return a list of dictionaries containing:
    index, name, alignment_filename, tree_filename
    """
    input_entries = []
    current_entry = {}
    
    with open(input_files, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if not line:  # Skip blank lines
            if current_entry:
                input_entries.append(current_entry)
                current_entry = {}
            continue
        
        if line.startswith('index:'):
            current_entry['index'] = line.split(':', 1)[1].strip()
        elif line.startswith('name:'):
            current_entry['name'] = line.split(':', 1)[1].strip()
        elif line.startswith('alignment_filename:'):
            current_entry['alignment_filename'] = line.split(':', 1)[1].strip()
        elif line.startswith('tree_filename:'):
            current_entry['tree_filename'] = line.split(':', 1)[1].strip()
    
    # Don't forget the last entry
    if current_entry:
        input_entries.append(current_entry)
    
    return input_entries

def parse_metrics_file(metrics_file):
    """
    Parse the metrics CSV file and return a list of dictionaries containing:
    Metric_name, phykit_command, metric_type
    """
    metrics = []
    
    with open(metrics_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Strip whitespace from keys and values
            metric = {
                'name': row['Metric_name'].strip(),
                'command': row[' phykit_command'].strip(),
                'type': row[' metric_type'].strip()
            }
            metrics.append(metric)
    
    return metrics

def expand_path(path_str):
    """Expand ~ and environment variables in path"""
    return os.path.expanduser(os.path.expandvars(path_str))

def run_phykit_metric(metric, alignment_file, tree_file):
    """
    Run a single phykit metric and return the raw output as a string.
    Returns error message if the metric cannot be run or an error occurs.
    """
    try:
        if metric['type'].lower() == 'alignment':
            cmd = ['phykit', metric['command'], alignment_file]
        elif metric['type'].lower() == 'tree':
            cmd = ['phykit', metric['command'], tree_file]
        elif metric['type'].lower() == 'both':
            cmd = ['phykit', metric['command'], '-a', alignment_file, '-t', tree_file]
        else:
            return "ERROR: Unknown metric type"
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0:
            output = result.stdout.strip()
            return output if output else "No output"
        else:
            return f"ERROR: {result.stderr.strip()}"
    except subprocess.TimeoutExpired:
        return "ERROR: Execution timeout"
    except Exception as e:
        return f"ERROR: {str(e)}"

def main():
    parser = argparse.ArgumentParser(
        description='run phykit metrics on input alignments and trees',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "input_files: text file with list of alignments and trees\n"
            "metrics: CSV file with each row as a metric to run. The file has three columns: Metric_name, phykit_command, and metric_type (tree, alignment, or both)\n"
            "output_filename: file to save output\n"
            "Example:\n"
            "  phykit_comparisons.py input_files.txt metrics.csv output.txt\n"
        )
    )
    parser.add_argument('input_files', help='text file with list of alignment and tree files to run phykit metrics on')
    parser.add_argument('metrics_file', help='CSV file with phykit metrics to compute')
    parser.add_argument('output_filename', help='file to save output')

    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output_filename), exist_ok=True)
    
    # Parse input files and metrics
    input_entries = parse_input_files(args.input_files)
    metrics = parse_metrics_file(args.metrics_file)
    
    # List to store results: [(metric_name, entry_name, output), ...]
    results = []
    
    # Run each metric for each input entry
    for entry in input_entries:
        alignment_file = expand_path(entry['alignment_filename'])
        tree_file = expand_path(entry['tree_filename'])
        
        # Check that files exist
        if not os.path.exists(alignment_file):
            print(f"Warning: Alignment file not found: {alignment_file}")
            continue
        if not os.path.exists(tree_file):
            print(f"Warning: Tree file not found: {tree_file}")
            continue
        
        print(f"Processing {entry['name']}...")
        
        for metric in metrics:
            print(f"  Running {metric['name']}...", end=' ')
            output = run_phykit_metric(metric, alignment_file, tree_file)
            results.append((metric['name'], entry['name'], output))
            print(f"Done")
    
    # Write output file
    output_file = args.output_filename
    with open(output_file, 'w') as f:
        f.write("PhyKIT Metrics Summary\n")
        f.write("=" * 80 + "\n\n")
        
        # Write metric descriptions
        f.write("Metrics Computed:\n")
        for i, metric in enumerate(metrics, 1):
            f.write(f"{i}. {metric['name']} ({metric['type']}): phykit {metric['command']}\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("Results\n")
        f.write("=" * 80 + "\n\n")
        
        # Sort results by metric name to group outputs by metric
        results_sorted = sorted(results, key=lambda x: x[0])
        
        # Write results with labeled headers, grouped by metric
        current_metric = None
        for metric_name, entry_name, output in results_sorted:
            if metric_name != current_metric:
                if current_metric is not None:
                    f.write("\n")
                f.write(f"Metric: {metric_name}\n")
                f.write("=" * 40 + "\n\n")
                current_metric = metric_name
            
            f.write(f"Input: {entry_name}\n")
            f.write("-" * 40 + "\n")
            f.write(output + "\n")
            f.write("\n")
    
    print(f"\nResults written to: {output_file}")

if __name__ == "__main__":
    main()