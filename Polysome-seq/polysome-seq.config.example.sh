# Copy to polysome-seq.config.sh and edit the paths.

# Input / output
# Put Input and Polysome paired-end FASTQ in the same directory:
#   huh7_input.R1.fastq.gz / huh7_input.R2.fastq.gz
#   huh7_poly.R1.fastq.gz  / huh7_poly.R2.fastq.gz
RAW_DIR="/path/to/polysome/raw_fastq"
OUTPUT_DIR="/path/to/polysome_output"

# References
STAR_INDEX="/path/to/STAR_genome_index"
GTF="/path/to/annotation.gtf"

# Threads and STAR memory
THREADS=8
STAR_SJDB_OVERHANG=149
STAR_MISMATCHES=2
STAR_BAM_SORT_RAM=64000000000

# PE RNA-seq-style length filter (applied by fastp to each mate)
MIN_LENGTH=20
MAX_LENGTH=200

# BAM filtering: keep reverse-strand RNA; only drop unmapped (flag 4)
MIN_MAPQ=10
FILTER_FLAGS=4

# featureCounts: reverse-stranded PE fragments, gene-level exon union
FC_FEATURE="exon"
FC_ATTRIBUTE="gene_id"
FC_STRANDNESS=2

# Yeast spike-in: host STAR unmapped PE pairs are remapped to yeast.
SPIKEIN=1
YEAST_STAR_INDEX="/path/to/yeast_STAR_genome_index"
YEAST_GTF="/path/to/yeast.gtf"
MIN_Y_READS=10000                 # warn if a library has fewer yeast Assigned fragments

# Manual TE. Sample sheet is optional if names contain input/poly keywords.
SAMPLE_SHEET=""                       # optional; copy polysome-seq.samples.example.tsv if names are not *input*/*poly*
INPUT_KEYWORD="input"
POLY_KEYWORD="poly"
TE_PSEUDOCOUNT=1
TE_MIN_TPM=1
# RSCRIPT="Rscript"
