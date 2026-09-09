# Copy to qti-seq.config.sh and edit the paths.

# Input / output
QTI_RAW_DIR="/path/to/qti/raw_fastq"     # paired-end: sample.R1.fastq.gz + sample.R2.fastq.gz
OUTPUT_DIR="/path/to/qti_output"
WITH_RNA=1                               # 1 = also process matched RNA-seq / Input libraries
RNA_RAW_DIR="/path/to/rna/raw_fastq"     # required when WITH_RNA=1

# References
SMRNA_INDEX="/path/to/rRNA_tRNA_snRNA_bowtie_index"   # Bowtie 1 prefix
STAR_INDEX="/path/to/STAR_genome_index"               # STAR genomeDir
GTF="/path/to/annotation.gtf"

# Threads and STAR memory
THREADS=8
STAR_SJDB_OVERHANG=149
STAR_MISMATCHES=2
STAR_BAM_SORT_RAM=64000000000

# QTI-seq length after PE merge (same RPF window as Ribo-seq)
MIN_LENGTH=18
MAX_LENGTH=40

# RNA-seq length on each mate (used only when WITH_RNA=1)
RNA_MIN_LENGTH=25
RNA_MAX_LENGTH=1000

# Bowtie rRNA/tRNA/snRNA depletion
BOWTIE_MISMATCHES=2
BOWTIE_MAX_ALIGNMENTS=1000

# BAM filtering
MIN_MAPQ=10
FILTER_FLAGS=1804   # unmapped + mate unmapped + secondary + QC fail + duplicate

# featureCounts (genome CDS counts; QTI libraries are single-end after merge)
FC_FEATURE="CDS"
FC_ATTRIBUTE="gene_id"
FC_STRANDNESS=0
