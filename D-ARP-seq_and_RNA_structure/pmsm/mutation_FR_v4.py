#!/usr/bin/env python3

# Compared with the previous version, this script also computes PMSM values.
# Usage:
#   python /home/pf/14T/scripts/mutation_FR_v4.py --input-dir bam_files --output-dir results

import os
import sys
import subprocess
import glob
import time
import datetime
from pathlib import Path
from multiprocessing import Pool, cpu_count
import argparse
import shutil

# Global variables
LOG_DIR = None
TIMESTAMP = None
MAIN_LOG = None
TIMING_LOG = None
PARALLEL_JOBS = 2
CLEANUP_TEMP = False
TRNA_FA = "genomic_UCSC.fa"
INPUT_DIR = None
OUTPUT_DIR = None
BASE_QUALITY = 0
MAPPING_QUALITY = 0


def init_logging(log_dir="./logs"):
    """Initialize the logging system."""
    global LOG_DIR, TIMESTAMP, MAIN_LOG, TIMING_LOG
    
    LOG_DIR = Path(log_dir)
    LOG_DIR.mkdir(exist_ok=True)
    
    TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    MAIN_LOG = LOG_DIR / f"rna_mutation_{TIMESTAMP}.log"
    TIMING_LOG = LOG_DIR / f"timing_{TIMESTAMP}.log"
    
    # Initialize the timing log
    with open(TIMING_LOG, 'w', encoding='utf-8') as f:
        f.write("File,Step,Duration,Status\n")


def log_message(message, file_log=None):
    """Write a log message."""
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    log_entry = f"[{timestamp}] {message}"
    
    # Print to the console
    print(log_entry)
    
    # Write to the main log
    with open(MAIN_LOG, 'a', encoding='utf-8') as f:
        f.write(log_entry + '\n')
    
    # Also write to a per-file log when provided
    if file_log:
        with open(file_log, 'a', encoding='utf-8') as f:
            f.write(log_entry + '\n')


def log_timing(entry):
    """Append a timing-log entry."""
    with open(TIMING_LOG, 'a', encoding='utf-8') as f:
        f.write(entry + '\n')
    print(entry)


def calculate_duration(start_time, end_time):
    """Format an elapsed duration as HH:MM:SS."""
    duration = int(end_time - start_time)
    hours = duration // 3600
    minutes = (duration % 3600) // 60
    seconds = duration % 60
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def complement(base):
    """Return the complementary base."""
    comp_map = {
        'A': 'T', 'a': 'T',
        'T': 'A', 't': 'A',
        'C': 'G', 'c': 'G',
        'G': 'C', 'g': 'C',
        'N': 'N', 'n': 'N'
    }
    return comp_map.get(base, base)


def process_mpileup_line(line):
    """Parse one mpileup line into strand-specific mutation counts."""
    parts = line.strip().split('\t')
    if len(parts) < 5:
        return []
    
    chr_name = parts[0]
    position = parts[1]
    ref_base = parts[2].upper()
    try:
        depth = int(parts[3])
    except (ValueError, IndexError):
        return []  # Skip invalid lines
    bases = parts[4] if depth > 0 else ""
    
    # Initialize plus- and minus-strand base counts
    base_types = ['A', 'C', 'G', 'T', 'N']
    count_plus = {b: 0 for b in base_types}
    count_minus = {b: 0 for b in base_types}
    plus_depth = 0
    minus_depth = 0
    
    # Initialize PMSM counts per strand
    P_plus, PM_plus, S_plus, SM_plus = 0, 0, 0, 0
    P_minus, PM_minus, S_minus, SM_minus = 0, 0, 0, 0
    
    # Parse the mpileup base string
    i = 0
    while i < len(bases):
        base = bases[i]
        
        # True when the next character is '$' (minus-strand stop)
        is_followed_by_end = (i + 1 < len(bases) and bases[i + 1] == '$')
        
        # Handle special mpileup characters
        if base == '^':
            # Read-start marker; only plus-strand starts are stop events
            i += 1  # Skip '^'
            if i < len(bases):
                # Skip the mapping-quality character
                i += 1
                if i < len(bases):
                    next_base = bases[i]
                    if next_base == '.':
                        # S: plus-strand read start matching the reference (RT stop)
                        S_plus += 1
                        count_plus[ref_base] += 1
                        plus_depth += 1
                        i += 1  # Skip the consumed base
                    elif next_base in 'ACGTN':
                        # SM: plus-strand read start with a mismatch (RT stop + mutation)
                        SM_plus += 1
                        count_plus[next_base] += 1
                        plus_depth += 1
                        i += 1  # Skip the consumed base
                    elif next_base == ',':
                        # Minus-strand read start; not a stop event
                        P_minus += 1
                        real_base = complement(ref_base)
                        count_minus[real_base] += 1
                        minus_depth += 1
                        i += 1  # Skip the consumed base
                    elif next_base in 'acgtn':
                        # Minus-strand read start with a mismatch; not a stop event
                        PM_minus += 1
                        real_base = complement(next_base)
                        count_minus[real_base] += 1
                        minus_depth += 1
                        i += 1  # Skip the consumed base
            continue
        elif base == '$':
            # Read-end marker; skip it (minus-strand stops are counted on the base)
            i += 1
            continue
        elif base in '+-':
            # Skip insertions/deletions
            i += 1
            len_str = ""
            while i < len(bases) and bases[i].isdigit():
                len_str += bases[i]
                i += 1
            if len_str:
                i += int(len_str)
            continue
        elif base == '*':
            # Deletion placeholder
            i += 1
            continue
        
        # Count bases, separating pass-through and stop events
        if base == '.':
            # Plus-strand match
            P_plus += 1
            count_plus[ref_base] += 1
            plus_depth += 1
        elif base == ',':
            if is_followed_by_end:
                # S: minus-strand read end matching the reference (RT stop)
                S_minus += 1
            else:
                # P: minus-strand pass-through match
                P_minus += 1
            # Minus-strand reads report the complement of the reference base
            real_base = complement(ref_base)
            count_minus[real_base] += 1
            minus_depth += 1
        elif base in 'ACGTN':
            # Plus-strand mismatch (not at a read start)
            PM_plus += 1
            count_plus[base] += 1
            plus_depth += 1
        elif base in 'acgtn':
            if is_followed_by_end:
                # SM: minus-strand read end with a mismatch (RT stop + mutation)
                SM_minus += 1
            else:
                # PM: minus-strand pass-through mismatch
                PM_minus += 1
            # mpileup reports plus-strand coordinates; convert to the sequenced base
            real_base = complement(base)
            count_minus[real_base] += 1
            minus_depth += 1
        
        i += 1
    
    results = []
    
    # Emit plus-strand counts when plus-strand reads cover this position
    if plus_depth > 0:
        row = [chr_name, position, ref_base, '+', str(plus_depth)]
        row.extend([str(count_plus[b]) for b in base_types])
        # Append PMSM values
        row.extend([str(P_plus), str(PM_plus), str(S_plus), str(SM_plus)])
        results.append('\t'.join(row))
    
    # Emit minus-strand counts when minus-strand reads cover this position
    # Base counts here are the actual sequenced bases
    if minus_depth > 0:
        row = [chr_name, position, ref_base, '-', str(minus_depth)]
        row.extend([str(count_minus[b]) for b in base_types])
        # Append PMSM values
        row.extend([str(P_minus), str(PM_minus), str(S_minus), str(SM_minus)])
        results.append('\t'.join(row))
    
    return results


def generate_mpileup(bam_file, prefix, file_log):
    """Generate an mpileup file from a BAM."""
    mpileup_file = OUTPUT_DIR / f"{prefix}-mpileup"
    log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Starting mpileup generation...", file_log)
    mpileup_start = time.time()
    
    try:
        # Run samtools mpileup
        cmd = [
            'samtools', 'mpileup',
            '-B',
            '-A',
            '-x',
            '-d', '100000000000',
            '-Q', str(BASE_QUALITY),  # Base-quality threshold
            '-q', str(MAPPING_QUALITY),  # Mapping-quality threshold
            '--ff', '0',  # Do not filter reads by FLAG (UNMAP, SECONDARY, QCFAIL, DUP)
            '-f', TRNA_FA,
            bam_file
        ]
        
        with open(str(mpileup_file), 'w', encoding='utf-8') as out_file:
            result = subprocess.run(cmd, stdout=out_file, stderr=subprocess.PIPE, text=True)
            
        if result.returncode != 0:
            raise Exception(f"Samtools mpileup failed: {result.stderr}")
        
        mpileup_end = time.time()
        mpileup_duration = calculate_duration(mpileup_start, mpileup_end)
        log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Mpileup completed. Duration: {mpileup_duration}", file_log)
        log_timing(f"{prefix},mpileup,{mpileup_duration},SUCCESS")
        
        # Log mpileup file size
        mpileup_size = os.path.getsize(str(mpileup_file)) / (1024 * 1024)  # MB
        log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Mpileup file size: {mpileup_size:.2f}MB", file_log)
        
        return mpileup_file, True
        
    except Exception as e:
        log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] ERROR: Mpileup generation failed! {str(e)}", file_log)
        log_timing(f"{prefix},mpileup,FAILED,ERROR")
        return None, False


def detect_mutations(mpileup_file, prefix, file_log):
    """Run strand-specific mutation calling on an mpileup file."""
    mutations_file = OUTPUT_DIR / f"{prefix}-mutations.txt"
    log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Starting strand-specific mutation detection...", file_log)
    log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Note: Negative strand will show actual sequenced bases", file_log)
    awk_start = time.time()
    
    try:
        # Create the output file and write the header
        with open(str(mutations_file), 'w', encoding='utf-8') as out_file:
            out_file.write("chr\tposition\tref_base\tstrand\tdepth\tA\tC\tG\tT\tN\tP\tPM\tS\tSM\n")
            
            # Process the mpileup file
            with open(str(mpileup_file), 'r', encoding='utf-8') as in_file:
                line_count = 0
                output_lines = 1
                plus_lines = 0
                minus_lines = 0
                for line in in_file:
                    line_count += 1
                    
                    # Process each line
                    results = process_mpileup_line(line)
                    for result in results:
                        out_file.write(result + '\n')
                        output_lines += 1
                        strand = result.split('\t', 4)[3]
                        plus_lines += strand == '+'
                        minus_lines += strand == '-'
                    
                    # Progress update every 10,000 lines
                    if line_count % 10000 == 0:
                        log_message(f"[Progress] Processed {line_count} lines...", file_log)
                
                log_message(f"[Progress] Finished processing {line_count} lines.", file_log)
                log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Input lines to process: {line_count}", file_log)
        
        awk_end = time.time()
        awk_duration = calculate_duration(awk_start, awk_end)
        log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Mutation detection completed. Duration: {awk_duration}", file_log)
        log_timing(f"{prefix},mutation_detection,{awk_duration},SUCCESS")
        
        # Log output-file stats
        output_size = os.path.getsize(str(mutations_file)) / (1024 * 1024)  # MB
        log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Output lines: {output_lines}", file_log)
        log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Output file size: {output_size:.2f}MB", file_log)
        
        log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Positive strand positions: {plus_lines}", file_log)
        log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Negative strand positions: {minus_lines}", file_log)
        
        return True
        
    except Exception as e:
        log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] ERROR: Mutation detection failed! {str(e)}", file_log)
        log_timing(f"{prefix},mutation_detection,FAILED,ERROR")
        return False


def detect_rna_mutations(bam_file):
    """Run strand-specific RNA mutation calling for one BAM file."""
    bam_path = Path(bam_file)
    prefix = bam_path.stem
    
    # Per-file log
    file_log = LOG_DIR / f"{prefix}_{TIMESTAMP}.log"
    overall_start = time.time()
    
    log_message("=" * 40, file_log)
    log_message(f"Processing: {prefix}", file_log)
    log_message(f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", file_log)
    log_message(f"Input file: {bam_file}", file_log)
    log_message(f"Reference: {TRNA_FA}", file_log)
    log_message("Analysis type: Single-end FR RNA-seq strand-specific mutation detection", file_log)
    log_message("Output format: Negative strand shows actual sequenced bases", file_log)
    log_message("=" * 40, file_log)
    
    # Step 1: generate mpileup
    mpileup_file, mpileup_success = generate_mpileup(bam_file, prefix, file_log)
    if not mpileup_success:
        return False
    
    # Step 2: strand-specific mutation calling
    detection_success = detect_mutations(mpileup_file, prefix, file_log)
    
    # Total runtime
    overall_end = time.time()
    total_duration = calculate_duration(overall_start, overall_end)
    
    log_message("=" * 40, file_log)
    log_message(f"File processing completed: {prefix}", file_log)
    log_message(f"Total processing time: {total_duration}", file_log)
    log_message("Output files:", file_log)
    log_message(f"  - Mutations: {prefix}-mutations.txt", file_log)
    log_message(f"End time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", file_log)
    log_message("=" * 40, file_log)
    
    log_timing(f"{prefix},TOTAL,{total_duration},COMPLETED")
    
    # Optional: remove intermediate files
    if CLEANUP_TEMP and mpileup_file:
        log_message(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] Cleaning up temporary files...", file_log)
        os.remove(str(mpileup_file))
    
    return detection_success


def get_system_info():
    """Collect basic system information."""
    info = {}
    info['cpu_cores'] = cpu_count()
    
    # Memory info
    try:
        result = subprocess.run(['free', '-h'], capture_output=True, text=True)
        for line in result.stdout.split('\n'):
            if line.startswith('Mem:'):
                info['memory'] = line.split()[1]
                break
    except:
        info['memory'] = 'Unknown'
    
    # Disk space
    try:
        result = subprocess.run(['df', '-h', '.'], capture_output=True, text=True)
        lines = result.stdout.split('\n')
        if len(lines) > 1:
            info['disk_space'] = lines[-1].split()[3]
    except:
        info['disk_space'] = 'Unknown'
    
    return info


def find_bam_files(input_dir, pattern_begin=None, pattern_mid=None):
    """Find BAM files matching optional name patterns.

    Args:
        input_dir (Path): Input directory
        pattern_begin (str): Filename prefix
        pattern_mid (str): Substring that must appear in the filename

    Returns:
        list: Matching BAM file paths
    """
    bam_files = []
    
    if pattern_begin:
        if pattern_mid:
            # Prefix match plus a required substring
            candidates = glob.glob(str(input_dir / f"{pattern_begin}*.bam"))
            bam_files = [f for f in candidates if pattern_mid in os.path.basename(f)]
        else:
            # Prefix match only
            bam_files = glob.glob(str(input_dir / f"{pattern_begin}*.bam"))
    else:
        if pattern_mid:
            # Substring match only
            candidates = glob.glob(str(input_dir / "*.bam"))
            bam_files = [f for f in candidates if pattern_mid in os.path.basename(f)]
        else:
            # All BAM files
            bam_files = glob.glob(str(input_dir / "*.bam"))
    
    return sorted(bam_files)


def check_requirements():
    """Check that the reference and samtools are available."""
    # Check the reference FASTA
    if not os.path.exists(TRNA_FA):
        log_message(f"ERROR: Reference file {TRNA_FA} not found!")
        return False
    log_message(f"Reference genome: {TRNA_FA}")
    
    # Check samtools
    if not shutil.which('samtools'):
        log_message("ERROR: samtools is not installed or not in PATH!")
        return False
    
    # Samtools version
    try:
        result = subprocess.run(['samtools', '--version'], capture_output=True, text=True)
        version = result.stdout.split('\n')[0]
        log_message(f"Samtools version: {version}")
    except:
        log_message("WARNING: Could not determine samtools version")
    
    return True


def create_results_readme():
    """Write RESULTS_README.txt describing the output format."""
    readme_content = """RNA Strand-Specific Mutation Detection Results
==============================================

Output Files Description:
-------------------------
*-mutations.txt
   - Complete mutation profile for each position
   - Format: chr position ref_base strand depth A C G T N P PM S SM
   - ref_base: reference base at this position (positive strand)
   - strand='+': mutations from positive strand genes (bases shown as in reference)
   - strand='-': mutations from negative strand genes (bases shown as actually sequenced)
   
   IMPORTANT FORMAT: 
   - Positive strand: A/C/G/T counts represent bases relative to reference genome
   - Negative strand: A/C/G/T counts represent ACTUAL SEQUENCED BASES (not reference frame)
   
   PMSM VALUES (RNA Modification Detection):
   - P: Pass - reads normally passing through this position (match reference)
   - PM: Pass & Mutation - reads normally passing through this position (show mutations)
   - S: Stop - RT stops at this position (positive strand: read start; negative strand: read end)
   - SM: Stop & Mutation - RT stops at this position with mutations
   
   RNA Modification Signal = (PM + SM + S) / Depth × 100%

Interpreting Results:
--------------------
- For positive strand genes (strand='+'):
  * Direct interpretation: counts show bases relative to reference
  * Example: ref=C, C=80, T=20 indicates potential C>T mutations
  
- For negative strand genes (strand='-'):
  * Direct interpretation: counts show actual sequenced bases from negative strand reads
  * Example: ref=G (complement=C in mRNA), C=80, T=20 shows actual C>T pattern in sequenced reads
  * No need for mental conversion - the counts show actual sequenced bases

RNA Modification Detection (PMSM Values):
-----------------------------------------
- P + PM = reads passing through this position (standard coverage)
- S + SM = reads where RT stops at this position (potential modification sites)
- Higher S/SM ratios indicate potential RNA modifications

STRAND-SPECIFIC STOP DETECTION:
- Positive strand: S/SM detected at read START positions (^quality.)/(^quality[ATCGN])
- Negative strand: S/SM detected at read END positions (,$)/([atcgn]$)
- This reflects the directionality of reverse transcription

Example Analysis:
  Position with depth=100, P=70, PM=10, S=15, SM=5
  * Passing reads: 80 (P+PM), mutation rate: 12.5% (PM/(P+PM))
  * Stopping reads: 20 (S+SM), mutation rate: 25% (SM/(S+SM))  
  * Potential modification signal: 30% ((PM+SM+S)/depth)

Quality Filters Applied:
------------------------
- Base quality: Configurable via -Q/--base-quality parameter (default: Q10)
- Mapping quality: Configurable via -q/--mapping-quality parameter (default: Q10)

Pipeline Steps:
---------------
1. Generate mpileup with quality filtering
2. Perform strand-specific mutation detection with base conversion

Note: 
- Positive strand coordinates and bases are in reference frame
- Negative strand coordinates are in reference frame, but bases show actual sequenced nucleotides
- This format makes it easier to directly observe mutation patterns in both strands
"""
    
    readme_file = OUTPUT_DIR / 'RESULTS_README.txt'
    with open(str(readme_file), 'w', encoding='utf-8') as f:
        f.write(readme_content)
    
    log_message(f"Results description saved to {readme_file}")


def generate_summary(batch_start, batch_end):
    """Print a batch processing summary."""
    batch_duration = calculate_duration(batch_start, batch_end)
    
    log_message("=" * 40)
    log_message("Processing Summary")
    log_message("=" * 40)
    log_message(f"Total processing time: {batch_duration}")
    log_message(f"Log files saved in: {LOG_DIR}")
    
    print("\n===== Processing Statistics =====")
    
    # Parse the timing log
    successful_files = []
    failed_files = []
    total_seconds = 0
    file_count = 0
    
    with open(TIMING_LOG, 'r', encoding='utf-8') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 4:
                if parts[3] == 'COMPLETED' and parts[1] == 'TOTAL':
                    successful_files.append(parts[0])
                    # Accumulate total seconds
                    time_parts = parts[2].split(':')
                    if len(time_parts) == 3:
                        total_seconds += int(time_parts[0]) * 3600 + int(time_parts[1]) * 60 + int(time_parts[2])
                        file_count += 1
                elif parts[3] == 'ERROR':
                    failed_files.append(parts[0])
    
    print("Successful files:")
    for f in sorted(set(successful_files)):
        print(f"  {f}")
    
    if failed_files:
        print("\nFailed files:")
        for f in sorted(set(failed_files)):
            print(f"  {f}")
    
    if file_count > 0:
        avg_seconds = total_seconds / file_count
        avg_time = calculate_duration(0, avg_seconds)
        print(f"\nAverage processing time per file:")
        print(f"  {avg_time} ({file_count} files processed)")
    
    log_message("=" * 40)
    log_message("All processing completed!")
    log_message(f"Main log: {MAIN_LOG}")
    log_message(f"Timing log: {TIMING_LOG}")
    log_message("=" * 40)


def process_all_bam_files(bam_files, parallel_jobs):
    """Process all BAM files, in parallel when requested."""
    log_message("Starting processing...")
    
    if parallel_jobs > 1 and len(bam_files) > 1:
        # Parallel processing
        with Pool(processes=min(parallel_jobs, len(bam_files))) as pool:
            pool.map(detect_rna_mutations, bam_files)
    else:
        # Serial processing
        for bam_file in bam_files:
            detect_rna_mutations(bam_file)


def main():
    """Command-line entry point."""
    global PARALLEL_JOBS, CLEANUP_TEMP, TRNA_FA, INPUT_DIR, OUTPUT_DIR, BASE_QUALITY, MAPPING_QUALITY
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='RNA Strand-Specific Mutation Detection Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage Examples:
  python mutation_FR_v3.py --input-dir bam_files --output-dir results
  python mutation_FR_v3.py --input-dir data --pattern-begin "sample" --parallel-jobs 4
  python mutation_FR_v3.py --input-dir data --pattern-mid "treated" --cleanup-temp
  python mutation_FR_v3.py --input-dir data --pattern-begin "ribo" --pattern-mid "control"
  python mutation_FR_v3.py --input-dir data -Q 20 -q 30 --output-dir results_high_quality
  python mutation_FR_v3.py --input-dir data --base-quality 15 --mapping-quality 20

File Filtering:
  --pattern-begin: Match files starting with specified prefix (e.g., "sample" matches "sample01.bam", "sample02.bam")
  --pattern-mid: Match files containing specified string (e.g., "treated" matches "sample_treated_1.bam", "data_treated.bam")
  Both patterns can be used together for more specific filtering
        """
    )
    parser.add_argument('--input-dir', type=str, default='.',
                        help='Input directory containing BAM files (default: current directory)')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory for results (default: current directory)')
    parser.add_argument('--parallel-jobs', type=int, default=2,
                        help='Number of parallel jobs (default: 2)')
    parser.add_argument('--log-dir', type=str, default=None,
                        help='Directory for log files (default: same as output directory)')
    parser.add_argument('--cleanup-temp', action='store_true',
                        help='Clean up temporary files after processing')
    parser.add_argument('--reference', type=str, default='genomic_UCSC.fa',
                        help='Reference genome file (default: genomic_UCSC.fa)')
    parser.add_argument('--pattern-begin', type=str, default=None,
                        help='File name prefix pattern for filtering specific samples')
    parser.add_argument('--pattern-mid', type=str, default=None,
                        help='File name middle pattern for filtering specific samples')
    parser.add_argument('-Q', '--base-quality', type=int, default=0,
                        help='Base quality threshold for mpileup (default: 0)')
    parser.add_argument('-q', '--mapping-quality', type=int, default=0,
                        help='Mapping quality threshold for mpileup (default: 0)')
    
    args = parser.parse_args()
    
    # Resolve input and output directories
    INPUT_DIR = Path(args.input_dir).resolve()
    OUTPUT_DIR = Path(args.output_dir).resolve()
    
    # Create the output directory
    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    
    # Set global options
    PARALLEL_JOBS = args.parallel_jobs
    CLEANUP_TEMP = args.cleanup_temp
    TRNA_FA = args.reference
    BASE_QUALITY = args.base_quality
    MAPPING_QUALITY = args.mapping_quality
    
    # Initialize logging; default to <output_dir>/logs when --log-dir is omitted
    log_dir = args.log_dir if args.log_dir else os.path.join(OUTPUT_DIR, "logs")
    init_logging(log_dir)
    
    # Start the main run
    log_message("=" * 40)
    log_message("RNA Strand-Specific Mutation Detection Pipeline")
    log_message("=" * 40)
    log_message("Analysis type: Single-end FR-oriented RNA-seq")
    log_message("Output format: Negative strand shows actual sequenced bases")
    log_message("Pipeline steps:")
    log_message("  1. Generate mpileup with quality filtering")
    log_message("  2. Perform strand-specific mutation detection")
    log_message("=" * 40)
    
    # Check runtime requirements
    if not check_requirements():
        sys.exit(1)
    
    # Collect BAM files
    bam_files = find_bam_files(INPUT_DIR, args.pattern_begin, args.pattern_mid)
    if not bam_files:
        if args.pattern_begin or args.pattern_mid:
            pattern_info = []
            if args.pattern_begin:
                pattern_info.append(f"prefix: {args.pattern_begin}")
            if args.pattern_mid:
                pattern_info.append(f"contains: {args.pattern_mid}")
            log_message(f"ERROR: No BAM files found matching patterns ({', '.join(pattern_info)}) in input directory: {INPUT_DIR}")
        else:
            log_message(f"ERROR: No BAM files found in input directory: {INPUT_DIR}")
        sys.exit(1)
    
    log_message(f"Input directory: {INPUT_DIR}")
    log_message(f"Output directory: {OUTPUT_DIR}")
    
    # Report filename filters
    if args.pattern_begin or args.pattern_mid:
        pattern_info = []
        if args.pattern_begin:
            pattern_info.append(f"prefix: {args.pattern_begin}")
        if args.pattern_mid:
            pattern_info.append(f"contains: {args.pattern_mid}")
        log_message(f"File filtering patterns: {', '.join(pattern_info)}")
    else:
        log_message("File filtering: processing all BAM files")
    
    log_message(f"Found {len(bam_files)} BAM files to process")
    log_message(f"Parallel jobs: {PARALLEL_JOBS}")
    
    # Log system information
    sys_info = get_system_info()
    log_message("System information:")
    log_message(f"  CPU cores: {sys_info['cpu_cores']}")
    log_message(f"  Memory: {sys_info.get('memory', 'Unknown')}")
    log_message(f"  Disk space: {sys_info.get('disk_space', 'Unknown')}")
    
    # Batch start time
    batch_start = time.time()
    
    # Process all BAM files
    process_all_bam_files(bam_files, PARALLEL_JOBS)
    
    # Batch end time
    batch_end = time.time()
    
    # Write the processing summary
    generate_summary(batch_start, batch_end)
    
    # Write RESULTS_README.txt
    create_results_readme()


if __name__ == "__main__":
    main()
