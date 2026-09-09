#!/usr/bin/env python3
"""
================================================================================
Script: merge_stop_mutation_v2.py
Purpose: Merge stop and mutation tables (pandas, performance-oriented)
================================================================================

Background:
  Join RT-stop bedGraph files from call_stop.py with mutation tables from
  mutation_FR, then compute mutation rate and RNA-modification Signal for
  downstream site annotation.

  Typical inputs:
    - call_stop.py -> *.bedgraph (stop signal)
    - mutation_FR_v4/v5.py -> *-mutations.txt (PMSM counts)

  Downstream uses:
    - {prefix}-merged.txt (stop_value + mutation + signal_value)
    - Input for annotation tools such as annotate_mRNA_site_ultimate_v3.py
    - Modification-site filtering by Signal

  Merge key: (chr, position, strand)
    - stop bedGraph end (0-based) = mutation position (1-based)
    - Outer merge: keep stop-only and mutation-only sites

  Derived metrics:
    - mutation_rate = (effective_depth - ref_count) / effective_depth
    - signal_value  = (PM + SM + S) / Depth x 100%

  v2 performance:
    - pandas reads plus vectorized math
    - bulk to_csv writes
    - --large-file-mode batches huge mutation files

Workflow:
  parse args -> init_logging -> read_fai_lengths
    -> find_paired_files (sample prefix, several bedGraph suffixes)
    -> per pair: merge_stop_mutation_optimized or process_large_file_in_batches
    -> finalize_merged_dataframe (fillna, mutation_rate, signal_value, chr_length)
    -> save {prefix}-merged.txt

PMSM columns:
    P  = Pass (match)            PM = Pass + mutation
    S  = Stop (RT stop, match)   SM = Stop + mutation

File pairing:
    stop:     {prefix}[_shifted_combined_sorted].bedgraph
    mutation: {prefix}-mutations.txt
    output:   {prefix}-merged.txt

CLI:
    --stop-dir              stop bedGraph directory (default: 8.stop)
    --mutation-dir          mutation directory (default: 9.mutation)
    --output-dir            output directory (default: 10.merged)
    --fai-file              chromosome-length file (.fai or .txt)
    --large-file-mode       batch oversized mutation files
    --batch-size            batch size (default: 1000000)

Notes:
  1. A later v3 (Polars) build is column-compatible; prefer it for huge tables
  2. Outer join keeps one-sided sites and fills defaults
  3. bedGraph end must correspond to mutation position in the same coordinate system
  4. Missing PMSM columns are filled with 0
  5. --large-file-mode batches mutations first, then adds stop-only sites

================================================================================
"""

import os
import pandas as pd
import numpy as np
import glob
from collections import defaultdict
import argparse
import warnings
import logging
import datetime
from pathlib import Path
warnings.filterwarnings('ignore')

# Global log-file path
_log_path = None


def init_logging(output_dir):
    """Initialize the logging system."""
    global _log_path
    
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(exist_ok=True, parents=True)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    _log_path = log_dir / f"merge_stop_mutation_{timestamp}.log"


def log_message(message):
    """Write a log message to the console and log file."""
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    log_entry = f"[{timestamp}] {message}"
    
    # Print to the console
    print(log_entry)
    
    # Write to the log file
    if _log_path:
        try:
            with open(_log_path, 'a', encoding='utf-8') as f:
                f.write(log_entry + '\n')
        except Exception as e:
            print(f"Warning: could not write log file: {e}")


def read_fai_lengths(fai_file):
    """
    Read chromosome lengths; infer the format from the file extension.
    - .fai: standard FASTA index, first two columns (name, length)
    - .txt: two-column table (name, length)
    - other extensions: detect column count
    """
    import os
    
    try:
        file_ext = os.path.splitext(fai_file)[1].lower()
        
        if file_ext == '.fai':
            # Standard .fai: at least 5 columns, use the first two
            log_message("  Detected .fai format; reading as a standard FASTA index")
            fai_df = pd.read_csv(fai_file, sep='\t', header=None, 
                                names=['seq_name', 'length', 'offset', 'linebases', 'linewidth'],
                                usecols=[0, 1])
        elif file_ext == '.txt':
            # Two-column .txt
            log_message("  Detected .txt format; reading two columns")
            fai_df = pd.read_csv(fai_file, sep='\t', header=None, 
                                names=['seq_name', 'length'])
        else:
            # Unknown extension: detect column count
            log_message(f"  Unknown extension {file_ext}; detecting column count...")
            with open(fai_file, 'r') as f:
                first_line = f.readline().strip()
                n_cols = len(first_line.split('\t'))
            
            if n_cols == 2:
                log_message(f"  Detected {n_cols} columns; reading two-column format")
                fai_df = pd.read_csv(fai_file, sep='\t', header=None, 
                                    names=['seq_name', 'length'])
            elif n_cols >= 5:
                log_message("Reading as a standard FASTA index")
                fai_df = pd.read_csv(fai_file, sep='\t', header=None, 
                                    names=['seq_name', 'length', 'offset', 'linebases', 'linewidth'],
                                    usecols=[0, 1])
            else:
                raise ValueError(f"Unsupported file format, column count: {n_cols}")
            
        return dict(zip(fai_df['seq_name'], fai_df['length']))
    except Exception as e:
        log_message(f"Warning: failed to read chromosome-length file: {e}")
        return {}


def find_paired_files(stop_dir, mutation_dir):
    """
    Pair stop and mutation files that share a sample prefix.
    """
    stop_files = glob.glob(os.path.join(stop_dir, "*.bedgraph"))
    mutation_files = glob.glob(os.path.join(mutation_dir, "*-mutations.txt"))

    def normalize_stop_prefix(file_path):
        basename = os.path.basename(file_path)
        known_suffixes = [
            "_shifted_combined_sorted.bedgraph",
            "_shifted_combined.bedgraph",
            ".bedgraph",
        ]
        for suffix in known_suffixes:
            if basename.endswith(suffix):
                return basename[:-len(suffix)]
        return os.path.splitext(basename)[0]

    def normalize_mutation_prefix(file_path):
        basename = os.path.basename(file_path)
        suffix = "-mutations.txt"
        if basename.endswith(suffix):
            return basename[:-len(suffix)]
        return os.path.splitext(basename)[0]
    
    # Sample prefixes, matching default names from call_stop.py and mutation_FR_v4.py
    stop_prefixes = {}
    for file in stop_files:
        prefix = normalize_stop_prefix(file)
        stop_prefixes[prefix] = file
    
    mutation_prefixes = {}
    for file in mutation_files:
        prefix = normalize_mutation_prefix(file)
        mutation_prefixes[prefix] = file
    
    # Collect matched pairs
    paired_files = []
    for prefix in stop_prefixes:
        if prefix in mutation_prefixes:
            paired_files.append((prefix, stop_prefixes[prefix], mutation_prefixes[prefix]))
    
    return paired_files


def read_bedgraph_file_optimized(file_path):
    """
    Read a stop bedGraph with pandas.
    """
    try:
        df = pd.read_csv(file_path, sep='\t', header=None,
                        names=['chr', 'start', 'end', 'stop_value', 'strand'],
                        dtype={'chr': str, 'start': int, 'end': int, 
                              'stop_value': float, 'strand': str})
        df['position_1base'] = df['end']
        return df
    except Exception as e:
        log_message(f"  Warning: bedGraph read failed, using fallback: {e}")
        # Fallback for empty files or malformed tables
        return pd.DataFrame(columns=['chr', 'start', 'end', 'stop_value', 'strand', 'position_1base'])


def read_mutation_file_optimized(file_path):
    """
    Read a mutation table with pandas.
    Uses flexible dtypes so slightly irregular files still load.
    """
    try:
        # Inspect the column count first
        with open(file_path, 'r') as f:
            header_line = f.readline()
            sample_line = f.readline()
            if sample_line:
                n_cols = len(sample_line.strip().split('\t'))
            else:
                n_cols = 14  # Default column count
        
        log_message(f"  Mutation file has {n_cols} columns")
        
        # Choose column names from the count; do not force dtypes
        if n_cols >= 14:
            col_names = ['chr', 'position_1base', 'ref_base', 'strand', 'depth', 
                        'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM']
        else:
            col_names = ['chr', 'position_1base', 'ref_base', 'strand', 'depth', 
                        'A', 'C', 'G', 'T', 'N']
        
        # Let pandas infer dtypes
        log_message("  Reading mutation file...")
        read_cols = min(n_cols, len(col_names))
        df = pd.read_csv(file_path, sep='\t', skiprows=1,
                        names=col_names[:read_cols],
                        usecols=list(range(read_cols)))
        
        log_message(f"  Read {len(df):,} mutation rows")
        
        # Add missing PMSM columns as zeros
        for col in ['P', 'PM', 'S', 'SM']:
            if col not in df.columns:
                df[col] = 0
        
        # Coerce numeric columns, tolerating bad values
        try:
            # Convert numeric columns
            numeric_cols = ['position_1base', 'depth', 'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM']
            for col in numeric_cols:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
        except Exception as convert_error:
            log_message(f"  Warning: dtype conversion issue: {convert_error}")
        
        return df
    except Exception as e:
        log_message(f"  ERROR: failed to read mutation file: {e}")
        log_message("  Returning an empty DataFrame; mutation fields will be defaults.")
        return pd.DataFrame(columns=['chr', 'position_1base', 'ref_base', 'strand', 'depth',
                                    'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM'])


def calculate_mutation_rate_vectorized(df):
    """
    Vectorized mutation-rate calculation.
    """
    # Vectorized ref_count via numpy.select
    conditions = [
        df['ref_base'] == 'A',
        df['ref_base'] == 'C',
        df['ref_base'] == 'G',
        df['ref_base'] == 'T'
    ]
    choices = [df['A'], df['C'], df['G'], df['T']]
    ref_count = np.select(conditions, choices, default=0)
    
    # Effective depth
    effective_depth = df['depth'] - df['N']
    
    # Mask of rows with usable counts
    total_bases = df['A'] + df['C'] + df['G'] + df['T']
    valid_mask = (df['depth'] > 0) & (effective_depth > 0) & (total_bases > 0)
    
    # Initialize mutation_rate
    df['mutation_rate'] = 0.0
    
    # Vectorized mutation rate
    df.loc[valid_mask, 'mutation_rate'] = (
        (effective_depth[valid_mask] - ref_count[valid_mask]) / effective_depth[valid_mask]
    )
    
    # Clip mutation rate to [0, 1]
    df['mutation_rate'] = df['mutation_rate'].clip(lower=0, upper=1)
    
    return df


def finalize_merged_dataframe(merged_df, chromosome_lengths=None):
    """
    Fill missing values and compute derived output columns.
    """
    # Bulk fill of missing values
    log_message("  Filling missing values...")
    # Dictionary fill is faster than per-column fillna
    fill_values = {
        'stop_value': 0,
        'ref_base': 'N',
        'depth': 0,
        'A': 0,
        'C': 0,
        'G': 0,
        'T': 0,
        'N': 0,
        'P': 0,
        'PM': 0,
        'S': 0,
        'SM': 0
    }
    for col, value in fill_values.items():
        if col not in merged_df.columns:
            merged_df[col] = value
    merged_df = merged_df.fillna(fill_values)

    if 'stop_value' in merged_df.columns:
        merged_df['stop_value'] = pd.to_numeric(merged_df['stop_value'], errors='coerce').fillna(0)
    
    # Cast integer columns
    int_cols = ['depth', 'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM', 'position']
    for col in int_cols:
        if col in merged_df.columns:
            merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce').fillna(0).astype(np.int32)
    
    log_message("  Computing mutation rate (vectorized)...")
    merged_df = calculate_mutation_rate_vectorized(merged_df)
    
    log_message("  Computing Signal (vectorized)...")
    # Signal = (PM + SM + S) / Depth x 100%
    signal_sum = merged_df['PM'] + merged_df['SM'] + merged_df['S']
    depth_mask = merged_df['depth'] > 0
    merged_df['signal_value'] = 0.0
    merged_df.loc[depth_mask, 'signal_value'] = (
        signal_sum[depth_mask] / merged_df.loc[depth_mask, 'depth'] * 100.0
    )
    
    # Chromosome length (vectorized map)
    if chromosome_lengths:
        log_message("  Adding chromosome lengths...")
        merged_df['chr_length'] = merged_df['chr'].map(chromosome_lengths).fillna(0).astype(np.int32)
    else:
        merged_df['chr_length'] = 0
    
    # Sort output rows
    log_message("  Sorting...")
    merged_df = merged_df.sort_values(['chr', 'position', 'strand'], ignore_index=True)
    
    return merged_df


def merge_stop_mutation_optimized(stop_df, mutation_df, chromosome_lengths=None):
    """
    Outer-merge stop and mutation tables.
    """
    # Rename position columns for the join
    stop_df_renamed = stop_df.rename(columns={'position_1base': 'position'})
    mutation_df_renamed = mutation_df.rename(columns={'position_1base': 'position'})
    
    # pandas outer join
    log_message("  Merging tables...")
    merged_df = pd.merge(
        mutation_df_renamed, 
        stop_df_renamed[['chr', 'position', 'strand', 'stop_value']], 
        on=['chr', 'position', 'strand'], 
        how='outer'
    )
    
    return finalize_merged_dataframe(merged_df, chromosome_lengths)


def save_merged_file_optimized(merged_df, output_path):
    """
    Write the merged table with pandas to_csv.
    """
    # Fixed output column order
    columns = ['chr', 'position', 'ref_base', 'strand', 'stop_value', 
               'depth', 'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM',
               'mutation_rate', 'signal_value', 'chr_length']
    
    # Round numeric columns
    output_df = merged_df[columns].copy()
    output_df['mutation_rate'] = output_df['mutation_rate'].round(4)
    output_df['signal_value'] = output_df['signal_value'].round(2)
    output_df['stop_value'] = output_df['stop_value'].astype(float)
    
    # Single write
    output_df.to_csv(output_path, sep='\t', index=False, float_format='%.4g')


def process_large_file_in_batches(stop_file, mutation_file, chromosome_lengths, batch_size=1000000):
    """
    Batch oversized mutation files to reduce peak memory.
    """
    # Stop file is usually small enough to load fully
    stop_df = read_bedgraph_file_optimized(stop_file)
    stop_df_renamed = stop_df.rename(columns={'position_1base': 'position'})
    stop_lookup = stop_df_renamed[['chr', 'position', 'strand', 'stop_value']]
    
    # Batch the mutation file
    merged_chunks = []
    mutation_key_chunks = []
    
    try:
        # Count rows for progress reporting
        with open(mutation_file, 'r') as f:
            header_line = f.readline()
            sample_line = f.readline()
            n_cols = len(sample_line.strip().split('\t')) if sample_line else 14
            total_lines = 1 + sum(1 for line in f) if sample_line else 0
        
        n_batches = (total_lines + batch_size - 1) // batch_size
        col_names = ['chr', 'position_1base', 'ref_base', 'strand', 'depth',
                     'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM']
        read_cols = min(n_cols, len(col_names))
        
        # Read and process in chunks
        for i, chunk in enumerate(pd.read_csv(mutation_file, sep='\t', skiprows=1, 
                                             chunksize=batch_size, header=None,
                                             names=col_names[:read_cols],
                                             usecols=list(range(read_cols)))):
            log_message(f"    Processing batch {i+1}/{n_batches}...")
            
            for col in ['P', 'PM', 'S', 'SM']:
                if col not in chunk.columns:
                    chunk[col] = 0
            
            chunk_renamed = chunk.rename(columns={'position_1base': 'position'})
            mutation_key_chunks.append(chunk_renamed[['chr', 'position', 'strand']].copy())

            # Left-merge this batch; add stop-only sites after all batches.
            merged_chunk = pd.merge(
                chunk_renamed,
                stop_lookup,
                on=['chr', 'position', 'strand'],
                how='left'
            )
            merged_chunk = finalize_merged_dataframe(merged_chunk, chromosome_lengths)
            merged_chunks.append(merged_chunk)
    
    except Exception as e:
        log_message(f"    Batch processing failed, falling back to the full merge: {e}")
        mutation_df = read_mutation_file_optimized(mutation_file)
        return merge_stop_mutation_optimized(stop_df, mutation_df, chromosome_lengths)

    if mutation_key_chunks:
        all_mutation_keys = pd.concat(mutation_key_chunks, ignore_index=True).drop_duplicates()
        stop_only = pd.merge(
            stop_lookup,
            all_mutation_keys,
            on=['chr', 'position', 'strand'],
            how='left',
            indicator=True
        )
        stop_only = stop_only[stop_only['_merge'] == 'left_only'].drop(columns=['_merge'])
    else:
        stop_only = stop_lookup

    if not stop_only.empty:
        log_message("    Adding stop-only sites...")
        merged_chunks.append(finalize_merged_dataframe(stop_only, chromosome_lengths))
    
    # Concatenate batches
    if merged_chunks:
        final_df = pd.concat(merged_chunks, ignore_index=True)
        final_df = final_df.sort_values(['chr', 'position', 'strand'], ignore_index=True)
        return final_df
    else:
        return pd.DataFrame()


def main():
    parser = argparse.ArgumentParser(description='Merge stop and mutation files (pandas)')
    parser.add_argument('--stop-dir', default='8.stop', help='Directory of stop bedGraph files')
    parser.add_argument('--mutation-dir', default='9.mutation', help='Directory of mutation files')
    parser.add_argument('--output-dir', default='10.merged', help='Output directory')
    parser.add_argument('--fai-file', required=False, help='FASTA index (.fai) or two-column length file')
    parser.add_argument('--large-file-mode', action='store_true', help='Batch oversized mutation files')
    parser.add_argument('--batch-size', type=int, default=1000000, help='Batch size in large-file mode')
    
    args = parser.parse_args()
    
    # Ensure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Initialize logging
    init_logging(args.output_dir)
    
    # Optional chromosome lengths
    chromosome_lengths = None
    if args.fai_file:
        try:
            chromosome_lengths = read_fai_lengths(args.fai_file)
            log_message(f"Loaded {len(chromosome_lengths)} chromosome lengths from the fai file")
        except Exception as e:
            log_message(f"Warning: could not read fai file {args.fai_file}: {str(e)}")
            log_message("Continuing without chromosome lengths")
    
    # Find matched file pairs
    paired_files = find_paired_files(args.stop_dir, args.mutation_dir)
    
    log_message(f"Found {len(paired_files)} matched file pair(s)")
    if args.large_file_mode:
        log_message(f"Large-file mode enabled, batch size: {args.batch_size}")
    log_message("")
    
    # Counters
    total_processed = 0
    total_errors = 0
    
    for idx, (prefix, stop_file, mutation_file) in enumerate(paired_files, 1):
        log_message(f"[{idx}/{len(paired_files)}] Processing: {prefix}")
        log_message(f"  Stop file: {os.path.basename(stop_file)}")
        log_message(f"  Mutation file: {os.path.basename(mutation_file)}")
        
        try:
            output_file = os.path.join(args.output_dir, f"{prefix}-merged.txt")
            if args.large_file_mode:
                # Reuse the v3 SQLite disk backend; older chunking still concatenated every chunk.
                from merge_stop_mutation_v3 import process_large_file
                process_large_file(
                    stop_file, mutation_file, chromosome_lengths, output_file,
                    args.batch_size
                )
                log_message(f"  Output: {output_file}")
                total_processed += 1
                log_message("")
                continue
            else:
                # Regular mode: load both tables fully
                stop_df = read_bedgraph_file_optimized(stop_file)
                mutation_df = read_mutation_file_optimized(mutation_file)
                
                log_message(f"  Stop rows: {len(stop_df):,}")
                log_message(f"  Mutation rows: {len(mutation_df):,}")
                
                # Merge
                merged_df = merge_stop_mutation_optimized(stop_df, mutation_df, chromosome_lengths)
            
            log_message(f"  Merged rows: {len(merged_df):,}")
            
            # Save
            save_merged_file_optimized(merged_df, output_file)
            log_message(f"  Output: {output_file}")
            
            total_processed += 1
            
        except Exception as e:
            log_message(f"  Error: {str(e)}")
            total_errors += 1
        
        log_message("")  # Blank line between samples
    
    log_message("=" * 50)
    log_message("Processing finished.")
    log_message(f"Succeeded: {total_processed} file(s)")
    if total_errors > 0:
        log_message(f"Failed: {total_errors} file(s)")
    log_message(f"Results written to {args.output_dir}")
    log_message(f"Log file: {_log_path}")


if __name__ == "__main__":
    main()
