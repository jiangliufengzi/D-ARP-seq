#!/usr/bin/env python3
"""
================================================================================
Script: call_stop.py
Purpose: Extract read-end positions from BAM files for stop-site analysis
================================================================================

Background:
  Used in RNA structure probing (DMS-seq, SHAPE-seq), ribosome profiling, and
  translatomics workflows to extract 5' or 3' ends of aligned reads.

  Core idea:
  - 5' or 3' ends mark reverse-transcriptase stops or RNA modification sites
  - A coordinate shift corrects adapter/primer offset
  - Plus- and minus-strand coverage are counted separately

  Typical inputs:
  - Sorted BAM files from Bowtie2/STAR (or similar)
  - Deduplicated BAMs (for example *.rmdup.bam)

  Typical uses:
  - Stop-site analysis and visualization
  - The truncation arm of mutation-truncation analysis
  - RNA modification calling (m1A, m5C, and similar RT-stop events)
  - Ribosome A/P/E-site analysis (Ribo-seq)

  Dependencies:
  - bedtools (bamtobed, shift, genomecov) for the keep-temp path
  - A two-column genome-size file (chrNameLength.txt)

================================================================================
Workflow (execution order):

    Input: BAM directory + genome-size file
      -> parse arguments, set up logging under output_dir/logs/
      -> scan BAM files (--pattern / --pattern-begin / --pattern-mid)
      -> process_bam_parallel()  (serial if workers=1, else ThreadPoolExecutor)
           process_single_bam():
             default: one BAM pass that accumulates strand-specific ends
             --keep-temp: BAM->BED -> shift -> plus/minus genomecov -> merge/sort
      -> output {sample}_shifted_combined_sorted.bedgraph

Core functions:
  1. setup_logger()           - log to file and console
  2. get_duration()           - format elapsed time as HH:MM:SS.mmm
  3. run_command()            - run a shell command
  4. process_bam_parallel()   - scan files and dispatch work
     -> process_single_bam()  - per-BAM processing
        or count_endpoints_direct() when intermediate files are discarded

Input:
  BAM: standard *.bam, preferably sorted
  Genome file (chrNameLength.txt): tab-separated chrom\\tlength

Output:
  *_shifted_combined_sorted.bedgraph
    five tab-separated columns: chrom, start, end, count, strand
  logs/bedtools_processing_*.log

CLI:
  -i, --input-dir     BAM input directory (default: current directory)
  -o, --output-dir    output directory (default: current directory)
  -g, --genome-file   genome chrNameLength file (required)
  -w, --workers       thread count (default: 1)
  -p, --pattern       BAM glob (default: *.bam)
  --pattern-begin     filename prefix filter
  --pattern-mid       filename substring filter
  -e, --end-type      5 or 3 (default: 5)
  --plus-shift        plus-strand shift (default: 0)
  --minus-shift       minus-strand shift (default: 0)
  -k, --keep-temp     keep intermediate BED/bedGraph files

Notes:
  1. bedtools must be on PATH when --keep-temp is used
  2. The genome file must use the same chromosome names as the BAM
  3. Logs are written under output_dir/logs/
  4. ThreadPoolExecutor is used because the work is I/O bound

================================================================================
"""

import glob
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
from pathlib import Path
import time
from datetime import datetime, timedelta
import logging
from collections import Counter
import pysam

# ---------- Logging ----------
def setup_logger(log_file):
    """Configure file and console logging."""
    logger = logging.getLogger('bedtools_processing')
    logger.setLevel(logging.INFO)
    
    # File handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    # Formatter with millisecond precision (3 digits)
    class MillisecondFormatter(logging.Formatter):
        def formatTime(self, record, datefmt=None):
            ct = self.converter(record.created)
            if datefmt:
                s = time.strftime(datefmt, ct)
            else:
                t = time.strftime('%Y-%m-%d %H:%M:%S', ct)
                s = '%s.%03d' % (t, record.msecs)
            return s
    
    formatter = MillisecondFormatter('[%(asctime)s] %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

# ---------- Helpers ----------
def get_duration(start_time, end_time):
    """Format elapsed time as HH:MM:SS.mmm."""
    duration = end_time - start_time
    hours = int(duration // 3600)
    minutes = int((duration % 3600) // 60)
    seconds = duration % 60  # Keep the fractional seconds
    return f"{hours:02d}:{minutes:02d}:{seconds:06.3f}"

def run_command(cmd, check=True):
    """Run a shell command."""
    try:
        result = subprocess.run(cmd, shell=True if isinstance(cmd, str) else False, 
                              capture_output=True, text=True, check=check)
        return result
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Command failed: {' '.join(cmd) if isinstance(cmd, list) else cmd}\n"
                          f"stderr: {e.stderr}")


def count_endpoints_direct(bam_file, genome_file, output_file, end_type,
                           plus_shift, minus_shift):
    """Count strand-specific ends in one BAM pass, without BED intermediates."""
    chrom_lengths = {}
    with open(genome_file) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                chrom_lengths[fields[0]] = int(fields[1])

    counts = Counter()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_end is None:
                continue
            chrom = bam.get_reference_name(read.reference_id)
            chrom_len = chrom_lengths.get(chrom)
            if chrom_len is None:
                continue
            strand = '-' if read.is_reverse else '+'
            shift = minus_shift if strand == '-' else plus_shift
            shifted_start = max(0, min(chrom_len, read.reference_start + shift))
            shifted_end = max(0, min(chrom_len, read.reference_end + shift))
            if shifted_start >= shifted_end:
                continue
            if end_type == '5':
                pos = shifted_end - 1 if strand == '-' else shifted_start
            else:
                pos = shifted_start if strand == '-' else shifted_end - 1
            counts[(chrom, pos, strand)] += 1

    with open(output_file, "w") as out:
        for (chrom, pos, strand), count in sorted(
            counts.items(), key=lambda item: (item[0][0], item[0][1], item[0][2] == '-')
        ):
            out.write(f"{chrom}\t{pos}\t{pos + 1}\t{count}\t{strand}\n")
    return len(counts)

# ---------- Core processing ----------
def process_single_bam(bam_file, genome_file, output_dir, logger, file_index, total_files, 
                       end_type='5', plus_shift=-1, minus_shift=1, keep_temp=False):
    """
    Process one BAM file.

    Args:
      bam_file: BAM path
      genome_file: genome-size file path
      output_dir: output directory
      logger: logger
      file_index: current file index
      total_files: total number of files
      end_type: '5' or '3' end to count
      plus_shift: plus-strand shift (default -1)
      minus_shift: minus-strand shift (default 1)
      keep_temp: keep intermediate files (default False)
    """
    try:
        bam_path = Path(bam_file)
        prefix = bam_path.stem  # Filename without extension
        output_path = Path(output_dir)
        
        # Per-file start time
        file_start_time = time.time()
        
        logger.info("")
        logger.info(f"========== Processing file [{file_index}/{total_files}]: {bam_path.name} ==========")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Settings: count {end_type}' end, plus shift={plus_shift}, minus shift={minus_shift}")
        logger.info(f"Intermediate files: {'keep' if keep_temp else 'delete'}")

        combined_sorted = str(output_path / f"{prefix}_shifted_combined_sorted.bedgraph")
        if not keep_temp:
            logger.info("Counting ends in a single BAM pass (no intermediate BED/bedGraph)...")
            site_count = count_endpoints_direct(
                bam_file, genome_file, combined_sorted, end_type,
                plus_shift, minus_shift
            )
            file_end_time = time.time()
            logger.info(f"Endpoint sites: {site_count:,}")
            logger.info(f"Final output: {combined_sorted}")
            logger.info(f"File runtime: {get_duration(file_start_time, file_end_time)}")
            return True, prefix
        
        # Intermediate file paths
        temp_files = []
        
        # Step 1: BAM to BED
        step_start = time.time()
        logger.info("Start: BAM to BED...")
        bed_file = str(output_path / f"{prefix}.bed")
        temp_files.append(bed_file)
        cmd = f"bedtools bamtobed -i {bam_file} > {bed_file}"
        run_command(cmd)
        step_end = time.time()
        logger.info(f"Done: BAM to BED (elapsed: {get_duration(step_start, step_end)})")
        
        # Step 2: Shift
        step_start = time.time()
        logger.info(f"Start: shift (plus:{plus_shift}, minus:{minus_shift})...")
        shifted_bed = str(output_path / f"{prefix}_shifted.bed")
        temp_files.append(shifted_bed)
        cmd = f"bedtools shift -m {minus_shift} -p {plus_shift} -i {bed_file} -g {genome_file} > {shifted_bed}"
        run_command(cmd)
        step_end = time.time()
        logger.info(f"Done: shift (elapsed: {get_duration(step_start, step_end)})")
        
        # genomecov end option from end_type
        end_param = f"-{end_type}"  # -5 or -3
        
        # Step 3: Plus strand coverage
        step_start = time.time()
        logger.info(f"Start: plus-strand coverage ({end_type}' end)...")
        plus_bedgraph = str(output_path / f"{prefix}_shifted_plus.bedgraph")
        temp_files.append(plus_bedgraph)
        cmd = f"bedtools genomecov -bg -strand + {end_param} -i {shifted_bed} -g {genome_file} > {plus_bedgraph}"
        run_command(cmd)
        step_end = time.time()
        logger.info(f"Done: plus-strand coverage (elapsed: {get_duration(step_start, step_end)})")
        
        # Step 4: Minus strand coverage
        step_start = time.time()
        logger.info(f"Start: minus-strand coverage ({end_type}' end)...")
        minus_bedgraph = str(output_path / f"{prefix}_shifted_minus.bedgraph")
        temp_files.append(minus_bedgraph)
        cmd = f"bedtools genomecov -bg -strand - {end_param} -i {shifted_bed} -g {genome_file} > {minus_bedgraph}"
        run_command(cmd)
        step_end = time.time()
        logger.info(f"Done: minus-strand coverage (elapsed: {get_duration(step_start, step_end)})")
        
        # Step 5: Merge plus and minus strand with strand information
        step_start = time.time()
        logger.info("Start: merge plus/minus strands and add strand...")
        
        combined_bedgraph = str(output_path / f"{prefix}_shifted_combined.bedgraph")
        temp_files.append(combined_bedgraph)
        combined_sorted = str(output_path / f"{prefix}_shifted_combined_sorted.bedgraph")
        
        # Write plus-strand rows
        with open(plus_bedgraph, 'r') as f_plus, \
             open(combined_bedgraph, 'w') as f_out:
            for line in f_plus:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    f_out.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t+\n")
        
        # Append minus-strand rows
        with open(minus_bedgraph, 'r') as f_minus, \
             open(combined_bedgraph, 'a') as f_out:
            for line in f_minus:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    f_out.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t-\n")
        
        # Sort
        cmd = f"LC_COLLATE=C sort -k1,1 -k2,2n {combined_bedgraph} > {combined_sorted}"
        run_command(cmd)
        
        step_end = time.time()
        logger.info(f"Done: merge plus/minus strands (elapsed: {get_duration(step_start, step_end)})")
        
        # Delete intermediates unless --keep-temp
        if not keep_temp:
            logger.info("Cleanup: deleting intermediate files...")
            for temp_file in temp_files:
                try:
                    if Path(temp_file).exists():
                        Path(temp_file).unlink()
                        logger.debug(f"  deleted: {temp_file}")
                except Exception as e:
                    logger.warning(f"  failed to delete {temp_file}: {e}")
            logger.info("Cleanup: intermediate files deleted")
        else:
            logger.info("Keep: intermediate files retained")
        
        # Per-file runtime
        file_end_time = time.time()
        file_total_duration = get_duration(file_start_time, file_end_time)
        
        logger.info(f"Finished {bam_path.name}")
        logger.info(f"Final output: {combined_sorted}")
        logger.info(f"File runtime: {file_total_duration}")
        
        return True, prefix
        
    except Exception as e:
        # Record runtime even on failure
        file_end_time = time.time()
        file_total_duration = get_duration(file_start_time, file_end_time)
        logger.error(f"[ERROR] Failed while processing {bam_file}: {e}")
        logger.info(f"Runtime (failed): {file_total_duration}")
        return False, None

def process_bam_parallel(input_dir, output_dir, genome_file, workers=1, pattern='*.bam',
                         pattern_begin=None, pattern_mid=None,
                         end_type='5', plus_shift=-1, minus_shift=1, keep_temp=False):
    """
    Process BAM files in serial or parallel.

    Args:
      input_dir: input directory
      output_dir: output directory
      genome_file: genome-size file (chrNameLength.txt)
      workers: thread count
      pattern: glob used when pattern_begin and pattern_mid are both None
      pattern_begin: filename prefix filter
      pattern_mid: filename substring filter
      end_type: '5' or '3' end to count
      plus_shift: plus-strand shift
      minus_shift: minus-strand shift
      keep_temp: keep intermediate files
    """
    # Ensure the output directory exists
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Log to output_dir/logs
    logs_dir = Path(output_dir) / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = str(logs_dir / f"bedtools_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    logger = setup_logger(log_file)
    
    # Sorting locale
    os.environ['LC_COLLATE'] = 'C'
    
    # Script start time
    script_start = time.time()
    logger.info("========== Script started ==========")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Genome file: {genome_file}")
    logger.info(f"Settings: count {end_type}' end, plus shift={plus_shift}, minus shift={minus_shift}")
    logger.info(f"Intermediate files: {'keep' if keep_temp else 'delete'}")
    
    # Find BAM files with flexible name filters
    input_path = Path(input_dir)
    
    # Prefer prefix/substring filters when set
    if pattern_begin or pattern_mid:
        if pattern_begin:
            if pattern_mid:
                # Prefix and substring
                all_bams = list(input_path.glob(f'{pattern_begin}*.bam'))
                bam_files = [f for f in all_bams if pattern_mid in f.name]
                logger.info(f"Match mode: prefix='{pattern_begin}', contains='{pattern_mid}'")
            else:
                # Prefix only
                bam_files = list(input_path.glob(f'{pattern_begin}*.bam'))
                logger.info(f"Match mode: prefix='{pattern_begin}'")
        else:
            # Substring only
            all_bams = list(input_path.glob('*.bam'))
            bam_files = [f for f in all_bams if pattern_mid in f.name]
            logger.info(f"Match mode: contains='{pattern_mid}'")
    else:
        # Fallback glob
        bam_files = list(input_path.glob(pattern))
        logger.info(f"Match mode: glob='{pattern}'")
    
    if not bam_files:
        if pattern_begin or pattern_mid:
            logger.error(f"No BAM files matched the filters under {input_dir}")
            logger.error(f"Filters: prefix='{pattern_begin}', contains='{pattern_mid}'")
        else:
            logger.error(f"No BAM files matching {pattern} under {input_dir}")
        return
    
    bam_count = len(bam_files)
    logger.info(f"Found {bam_count} BAM file(s) to process")
    
    # List matched files (first 10)
    if bam_count > 0:
        logger.info("Matched files:")
        for i, bam_file in enumerate(bam_files[:10], 1):
            logger.info(f"  {i}. {bam_file.name}")
        if bam_count > 10:
            logger.info(f"  ... {bam_count - 10} more file(s) not shown")
    
    if workers == 1:
        # Serial processing
        for idx, bam_file in enumerate(bam_files, 1):
            process_single_bam(str(bam_file), genome_file, output_dir, logger, idx, bam_count,
                             end_type, plus_shift, minus_shift, keep_temp)
    else:
        # Parallel processing
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = []
            for idx, bam_file in enumerate(bam_files, 1):
                futures.append(
                    executor.submit(
                        process_single_bam,
                        str(bam_file),
                        genome_file,
                        output_dir,
                        logger,
                        idx,
                        bam_count,
                        end_type,
                        plus_shift,
                        minus_shift,
                        keep_temp
                    )
                )
            
            # Wait for all tasks
            success_count = 0
            fail_count = 0
            for future in as_completed(futures):
                try:
                    success, prefix = future.result()
                    if success:
                        success_count += 1
                    else:
                        fail_count += 1
                except Exception as e:
                    logger.error(f"[EXCEPTION] Sample processing failed: {e}")
                    fail_count += 1
            
            logger.info(f"Finished - success: {success_count}, failed: {fail_count}")
    
    # Script end time
    script_end = time.time()
    total_duration = get_duration(script_start, script_end)
    logger.info("")
    logger.info("========== Script finished ==========")
    logger.info(f"Total runtime: {total_duration}")
    logger.info(f"Log file: {log_file}")

# ---------- CLI entry point ----------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Process BAM files in parallel to call stop sites",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all BAM files in the current directory
  python call_stop.py

  # Set input, output, and genome-size file
  python call_stop.py -i /path/to/bam/files -o /path/to/output -g /path/to/genome.txt

  # Filter by filename prefix/substring
  python call_stop.py -o ./results --pattern-begin "sample" --pattern-mid "treated"
  python call_stop.py -o ./output --pattern-begin "ctrl"
  python call_stop.py -o ./processed --pattern-mid "rep1"

  # Glob pattern (used when --pattern-begin and --pattern-mid are omitted)
  python call_stop.py -o ./results --pattern "*.sorted.bam"

  # Parallel processing
  python call_stop.py -w 4 -o ./parallel_output --pattern-begin "sample"

  # Custom end type and shifts; keep intermediates
  python call_stop.py -o ./custom_output --end-type 3 --plus-shift 0 --minus-shift 0 --keep-temp
        """)
    parser.add_argument('--input-dir', '-i', default='.', help='BAM directory (default: current directory)')
    parser.add_argument('--output-dir', '-o', default='.', help='Output directory (default: current directory)')
    parser.add_argument('--genome-file', '-g', 
                       default='/home/pf/14T/index/human/hg38/human_hg38_genome_star/chrNameLength.txt',
                       help='Path to the genome chrNameLength file')
    parser.add_argument('--workers', '-w', type=int, default=1, 
                       help='Number of worker threads (default: 1, serial)')
    parser.add_argument('--pattern', '-p', default='*.bam',
                       help='BAM glob (default: *.bam; used when --pattern-begin and --pattern-mid are omitted)')
    parser.add_argument('--pattern-begin', type=str,
                       help='Filename prefix used to select samples')
    parser.add_argument('--pattern-mid', type=str,
                       help='Filename substring used to select samples')
    parser.add_argument('--end-type', '-e', choices=['5', '3'], default='5',
                       help="Which read end to count: '5' or '3' (default: 5)")
    parser.add_argument('--plus-shift', type=int, default=0,
                       help='Plus-strand shift (default: 0)')
    parser.add_argument('--minus-shift', type=int, default=0,
                       help='Minus-strand shift (default: 0)')
    parser.add_argument('--keep-temp', '-k', action='store_true',
                       help='Keep intermediate files (default: delete them and keep only the final output)')
    
    args = parser.parse_args()
    
    # Require the genome-size file
    if not Path(args.genome_file).exists():
        print(f"[ERROR] Genome file not found: {args.genome_file}")
        exit(1)
    
    # Run
    process_bam_parallel(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        genome_file=args.genome_file,
        workers=args.workers,
        pattern=args.pattern,
        pattern_begin=args.pattern_begin,
        pattern_mid=args.pattern_mid,
        end_type=args.end_type,
        plus_shift=args.plus_shift,
        minus_shift=args.minus_shift,
        keep_temp=args.keep_temp
    )
