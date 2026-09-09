#!/usr/bin/env bash
# Copy this file to config.sh and edit the paths for your dataset.

# Input FASTQ root. It is searched recursively. Files must be named
# SAMPLE.R1.fastq.gz and SAMPLE.R2.fastq.gz, with SAMPLE matching sample_id.
RAW_DIR="/path/to/1.raw_fq"

# Output directory. It is created automatically.
OUTPUT_DIR="/path/to/output"

# Bowtie 1 index prefix (do not include .1.ebwt). If it is absent, the
# pipeline builds it automatically from REFERENCE_FASTA.
BOWTIE_INDEX="/path/to/human_hg38_smrna_bowtie"

# The FASTA used to build the Bowtie index. Its .fai file is generated
# automatically with samtools faidx when absent.
REFERENCE_FASTA="/path/to/human_hg38_smrna.fa"

# Optional but strongly recommended. Set this to the SHA256 of the exact
# reference FASTA; the published HeLa workflow used:
# b0571c3a903633695f4a27a48d59ecf0b528852389f6e4fee91571e79db731d3
EXPECTED_REFERENCE_SHA256=""

# Experimental design for background correction and replicate merging.
# Start from samples.example.csv and keep sample_id identical to the FASTQ
# prefix (the part before .R1.fastq.gz).
SAMPLE_SHEET="/path/to/samples.csv"

# Resource settings.
THREADS=8
JAVA_MEMORY="16G"

# Read processing settings matching the current DHU-seq workflow.
FIRST_PASS_MIN_LENGTH=18
MIN_LENGTH=20
MAX_LENGTH=200
UMI_LOCATION="read2"
UMI_LENGTH=10
UMI_PREFIX="UMI"

# Bowtie and BAM filtering settings.
BOWTIE_MISMATCHES=3
BOWTIE_MAX_ALIGNMENTS=3
MIN_MAPQ=10
FILTER_FLAGS=1804

# UMICollapse settings.
UMI_EDIT_DISTANCE=1
UMI_ALGORITHM="adj"
UMI_SEPARATOR=":UMI_"
