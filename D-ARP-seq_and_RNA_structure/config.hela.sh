#!/usr/bin/env bash
# Ready-to-run configuration for reproducing the published HeLa folder 14.
# Put raw FASTQs under raw_fastq/ and the exact reference at the path below.

CONFIG_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

RAW_DIR="${CONFIG_DIR}/raw_fastq"
OUTPUT_DIR="${CONFIG_DIR}/results"
REFERENCE_FASTA="${CONFIG_DIR}/reference/human_hg38_smrna.fa"
BOWTIE_INDEX="${CONFIG_DIR}/reference/human_hg38_smrna_bowtie"
SAMPLE_SHEET="${CONFIG_DIR}/samples.hela.csv"

# Exact reference used by the existing HeLa folder-14 analysis.
EXPECTED_REFERENCE_SHA256="b0571c3a903633695f4a27a48d59ecf0b528852389f6e4fee91571e79db731d3"

THREADS=8
JAVA_MEMORY="16G"

FIRST_PASS_MIN_LENGTH=18
MIN_LENGTH=20
MAX_LENGTH=200
UMI_LOCATION="read2"
UMI_LENGTH=10
UMI_PREFIX="UMI"

BOWTIE_MISMATCHES=3
BOWTIE_MAX_ALIGNMENTS=3
MIN_MAPQ=10
FILTER_FLAGS=1804

UMI_EDIT_DISTANCE=1
UMI_ALGORITHM="adj"
UMI_SEPARATOR=":UMI_"
