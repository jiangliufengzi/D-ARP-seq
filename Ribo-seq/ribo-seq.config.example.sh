# Copy to ribo-seq.config.sh and edit the paths.

# Input / output
RIBO_RAW_DIR="/path/to/ribo/raw_fastq"   # paired-end: sample.R1.fastq.gz + sample.R2.fastq.gz
OUTPUT_DIR="/path/to/ribo_output"
WITH_RNA=1                               # 1 = also process matched RNA-seq / Input libraries
RNA_RAW_DIR="/path/to/rna/raw_fastq"     # required when WITH_RNA=1

# References
SMRNA_INDEX="/path/to/rRNA_tRNA_snRNA_bowtie_index"   # Bowtie 1 prefix
TX_INDEX="/path/to/longest_transcriptome_hisat2"      # HISAT2 prefix

# Threads
THREADS=8

# Ribo-seq length after PE merge (typical RPF window)
MIN_LENGTH=18
MAX_LENGTH=40

# RNA-seq length on each mate (used only when WITH_RNA=1)
RNA_MIN_LENGTH=25
RNA_MAX_LENGTH=1000

# Bowtie rRNA/tRNA/snRNA depletion
BOWTIE_MISMATCHES=2
BOWTIE_MAX_ALIGNMENTS=1000

# HISAT2 strandness. Single-end Ribo after merge uses F (from RF libraries).
# Paired RNA uses RF. Leave empty to omit --rna-strandness.
HISAT2_STRANDNESS="F"
HISAT2_RNA_STRANDNESS="RF"

# BAM filtering
MIN_MAPQ=10
FILTER_FLAGS=1804   # unmapped + mate unmapped + secondary + QC fail + duplicate

# Yeast spike-in (host-unmapped Ribo reads -> yeast transcriptome).
# RNA/Input is mapped to yeast from clean PE FASTQ, matching the original pipeline.
SPIKEIN=1
YEAST_TX_INDEX="/path/to/yeast_longest_hisat2"
YEAST_HISAT2_STRANDNESS="F"       # single-end Ribo after merge
YEAST_HISAT2_RNA_STRANDNESS="RF"  # paired RNA
MIN_Y_READS=1000                  # warn if a library has fewer yeast primary-mapped reads

# Manual TE (requires WITH_RNA=1). Copy ribo-seq.samples.example.tsv.
SAMPLE_SHEET="/path/to/samples.tsv"   # columns: sample, group, assay, pair
TE_PSEUDOCOUNT=1
TE_MIN_TPM=1
# RSCRIPT="Rscript"                  # optional; default is Rscript on PATH
