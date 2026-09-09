# Reproducible DHU-seq upstream pipeline

This repository contains a minimal workflow for processing paired-end DHU-seq
reads through background correction and biological-replicate merging. Standard
preprocessing steps are written directly in Bash. Existing PMSM, Sprinzl,
background-correction and replicate-merging implementations are included in
`pmsm/` and `postprocess/`.

## Workflow

1. `fastp`: paired-end QC, UMI extraction and read merging
2. `seqkit`: length filtering and reverse complementation
3. `bowtie`: alignment to the tRNA reference
4. `samtools` and `awk`: BAM sorting, indexing, MAPQ/FLAG filtering and
   removal of alignments containing insertions or deletions
5. `UMICollapse`: UMI-aware deduplication
6. Original `call_stop.py`, `mutation_FR_v4.py` and
   `merge_stop_mutation_v2.py`: PMSM signal generation
7. Original `annotate_tRNA_position_sprinzl.py`: Sprinzl annotation
8. Original `01_dhu_analysis.R`: background correction and DHU-site filtering
9. Original `02_dhu_analysis_merge_rep.R`: biological-replicate merging

## Input requirements

The FASTQ root is searched recursively. FASTQ files must use matching names:

```text
SAMPLE.R1.fastq.gz
SAMPLE.R2.fastq.gz
```

The user must provide:

- the exact reference FASTA (the Bowtie 1 index is built automatically);
- paired-end FASTQ files.
- a sample sheet defining background samples and replicate pairs.

The reference FASTA and Bowtie index must describe the same sequences. The
pipeline creates both the Bowtie 1 index and FASTA `.fai` automatically when
they are absent. `config.hela.sh` also checks the reference SHA256 before any
analysis starts.

## Installation

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
conda activate dhu-pmsm
```

## Exact HeLa reproduction

Place the inputs as follows:

```text
raw_fastq/SAMPLE.R1.fastq.gz
raw_fastq/SAMPLE.R2.fastq.gz
reference/human_hg38_smrna.fa
```

The FASTQ files may be in subdirectories of `raw_fastq/`. All 30 required
sample names are fixed in `samples.hela.csv`. The exact reference FASTA must
have this SHA256:

```text
b0571c3a903633695f4a27a48d59ecf0b528852389f6e4fee91571e79db731d3
```

Then the complete analysis is one command:

```bash
bash run_pipeline.sh config.hela.sh
```

The repository intentionally does not contain raw sequencing data or the
5-MB custom small-RNA reference FASTA. They must be deposited alongside the
code (for example in a GitHub Release, Zenodo or GEO/SRA-associated archive).
Code alone cannot reproduce folder 14 without these two input assets.

## Configuration for another experiment

```bash
cp config.example.sh config.sh
```

For a different experiment, start from the minimal example:

```bash
cp samples.example.csv samples.csv
```

Edit `config.sh`, especially `RAW_DIR`, `OUTPUT_DIR`, `BOWTIE_INDEX`,
`REFERENCE_FASTA` and `SAMPLE_SHEET`. In `samples.csv`, `sample_id` must be
identical to the corresponding FASTQ prefix before `.R1.fastq.gz`.

Sample-sheet columns:

- `sample_id`: unique FASTQ/sample prefix
- `role`: `background` or `treatment`
- `background_set`: connects each treatment to one or more backgrounds
- `merge_group`: output group shared by treatment replicate 1 and replicate 2
- `replicate`: replicate number (`1` or `2`)
- `summary_name`: prefix of the replicate-overlap summary file
- `summary_group`: condition label written into that summary

## Run

```bash
bash run_pipeline.sh config.sh
```

Final background-corrected, replicate-merged tables are written to
`OUTPUT_DIR/14_merge_DHU_site_rep/`.

Before publishing, verify that the copied original analysis code and Sprinzl
map have not changed:

```bash
sha256sum -c SOURCE_FILES.sha256
```

To compare a reproduced HeLa folder 14 with the existing four reference
folders, run:

```bash
python compare_14_results.py \
  --reference /path/to/HeLa_NaBH4/14_merge_DHU_site_rep \
  --reference /path/to/HeLa_ARP/14_merge_DHU_site_rep \
  --reference /path/to/HeLa_IP/14_merge_DHU_site_rep \
  --reference /path/to/HeLa_test/14_merge_DHU_site_rep \
  --reproduced results/14_merge_DHU_site_rep
```

## Output directories

```text
2.rmumi_fq/       first-pass merged/reverse-complemented reads
3.clean_fq/       second-pass filtered reads
4.mapped/         coordinate-sorted BAM files
5.filtered_bam/   MAPQ/flag/indel-filtered BAM files
6.bam_rmdup/      UMI-deduplicated BAM files
8.stop/           RT-stop signal
9.mutation/       mutation and PMSM components
10.merge/         merged PMSM tables
11.anno_sprinzl/  Sprinzl-annotated PMSM tables
13_calu_bg_tRNA/  background-corrected per-sample tables
14_merge_DHU_site_rep/ replicate-merged final tables
logs/             fastp and Bowtie logs
```

Folder 12 is intentionally omitted because the existing folder-14 workflow
uses all-tRNA results from folder 13, not the oligo-specific folder-12 branch.

## Important parameters

The HeLa configuration reproduces the current upstream settings: UMI in
read 2 with length 10, first-pass merged reads of 18–200 nt, second-pass reads
of 20–200 nt, Bowtie `-v 3 -m 3 --best
--strata --norc`, BAM filtering with MAPQ 10 and flag mask 1804, and
UMICollapse adjacency deduplication.

The post-processing runner also preserves the original HeLa-specific rules:
background correction, canonical D-loop positions, mitochondrial exceptions,
oligo removal, noncanonical-site filtering, two-replicate intersection, the
mitochondrial OR exception, Met-tRNA helper output and replicate-overlap
summary tables.
