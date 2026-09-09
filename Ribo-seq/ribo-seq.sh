#!/usr/bin/env bash
# Ribo-seq upstream workflow: raw paired-end FASTQ to transcriptome count matrices
# plus yeast spike-in size factors.
# QC (fastp merge) -> rRNA/tRNA depletion (Bowtie) -> HISAT2 transcriptome
# -> MAPQ/flag/length filter (samtools) -> transcript counts (samtools idxstats)
# -> host-unmapped Ribo (and clean RNA) to yeast HISAT2 -> size_factor = Y / geomean(Y)
# -> manual TE: classic TPM, then Ribo TPM x te_ratio(Q); RNA Q=1.
# Optional matched RNA-seq/Input libraries are kept paired-end.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/ribo-seq.config.sh}"

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: configuration file not found: $CONFIG" >&2
  echo "Copy ribo-seq.config.example.sh to ribo-seq.config.sh and edit the paths first." >&2
  exit 1
fi

# shellcheck source=/dev/null
source "$CONFIG"

WITH_RNA="${WITH_RNA:-0}"
SPIKEIN="${SPIKEIN:-0}"
HISAT2_STRANDNESS="${HISAT2_STRANDNESS:-}"
HISAT2_RNA_STRANDNESS="${HISAT2_RNA_STRANDNESS:-}"
YEAST_HISAT2_STRANDNESS="${YEAST_HISAT2_STRANDNESS:-}"
YEAST_HISAT2_RNA_STRANDNESS="${YEAST_HISAT2_RNA_STRANDNESS:-}"
MIN_Y_READS="${MIN_Y_READS:-1000}"
TE_PSEUDOCOUNT="${TE_PSEUDOCOUNT:-1}"
TE_MIN_TPM="${TE_MIN_TPM:-1}"
RSCRIPT="${RSCRIPT:-Rscript}"
SAMPLE_SHEET="${SAMPLE_SHEET:-}"

required_variables=(
  RIBO_RAW_DIR OUTPUT_DIR SMRNA_INDEX TX_INDEX THREADS
  MIN_LENGTH MAX_LENGTH BOWTIE_MISMATCHES BOWTIE_MAX_ALIGNMENTS
  MIN_MAPQ FILTER_FLAGS
)
if [[ "$WITH_RNA" == "1" ]]; then
  required_variables+=(RNA_RAW_DIR RNA_MIN_LENGTH RNA_MAX_LENGTH)
fi
if [[ "$SPIKEIN" == "1" ]]; then
  required_variables+=(YEAST_TX_INDEX)
fi
for variable in "${required_variables[@]}"; do
  if [[ -z "${!variable:-}" ]]; then
    echo "ERROR: required configuration variable is empty: $variable" >&2
    exit 1
  fi
done

for program in fastp seqkit bowtie hisat2 samtools gzip; do
  if ! command -v "$program" >/dev/null 2>&1; then
    echo "ERROR: required program is not available in PATH: $program" >&2
    exit 1
  fi
done

if [[ ! -d "$RIBO_RAW_DIR" ]]; then
  echo "ERROR: Ribo-seq FASTQ directory not found: $RIBO_RAW_DIR" >&2
  exit 1
fi
if [[ "$WITH_RNA" == "1" && ! -d "$RNA_RAW_DIR" ]]; then
  echo "ERROR: RNA-seq FASTQ directory not found: $RNA_RAW_DIR" >&2
  exit 1
fi
if [[ ! -f "${SMRNA_INDEX}.1.ebwt" && ! -f "${SMRNA_INDEX}.1.ebwtl" ]]; then
  echo "ERROR: Bowtie rRNA/tRNA index not found: ${SMRNA_INDEX}.1.ebwt[ l ]" >&2
  exit 1
fi
if [[ ! -f "${TX_INDEX}.1.ht2" && ! -f "${TX_INDEX}.1.ht2l" ]]; then
  echo "ERROR: HISAT2 transcriptome index not found: ${TX_INDEX}.1.ht2[ l ]" >&2
  exit 1
fi
if [[ "$SPIKEIN" == "1" && ! -f "${YEAST_TX_INDEX}.1.ht2" && ! -f "${YEAST_TX_INDEX}.1.ht2l" ]]; then
  echo "ERROR: yeast HISAT2 transcriptome index not found: ${YEAST_TX_INDEX}.1.ht2[ l ]" >&2
  exit 1
fi

sample_name_from_r1() {
  local base
  base="$(basename -- "$1")"
  base="${base%.R1.fastq.gz}"
  base="${base%.R1.fq.gz}"
  printf '%s\n' "$base"
}

find_r2() {
  local r1="$1"
  local dir sample
  dir="$(dirname -- "$r1")"
  sample="$(sample_name_from_r1 "$r1")"
  local candidate
  for candidate in \
    "${dir}/${sample}.R2.fastq.gz" \
    "${dir}/${sample}.R2.fq.gz"
  do
    if [[ -f "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

collect_r1_files() {
  local raw_dir="$1"
  local -a found=()
  local r1
  shopt -s nullglob
  for r1 in "${raw_dir}"/*.R1.fastq.gz "${raw_dir}"/*.R1.fq.gz; do
    [[ -e "$r1" ]] || continue
    found+=("$r1")
  done
  shopt -u nullglob
  if (( ${#found[@]} == 0 )); then
    echo "ERROR: no *.R1.fastq.gz or *.R1.fq.gz files found in $raw_dir" >&2
    exit 1
  fi
  printf '%s\n' "${found[@]}"
}

filter_sorted_bam() {
  local input_bam="$1"
  local output_bam="$2"
  local min_len="${3:-}"
  local max_len="${4:-}"
  samtools view \
    -@ "$THREADS" \
    -h \
    -q "$MIN_MAPQ" \
    -F "$FILTER_FLAGS" \
    "$input_bam" |
  awk -v min_len="$min_len" -v max_len="$max_len" '
    BEGIN { OFS="\t" }
    /^@/ { print; next }
    $6 ~ /[ID]/ { next }
    {
      n = length($10)
      if (min_len != "" && n < min_len + 0) next
      if (max_len != "" && n > max_len + 0) next
      print
    }
  ' |
  samtools view -@ "$THREADS" -b -o "$output_bam" -
  samtools index -@ "$THREADS" "$output_bam"
}

merge_idxstats_matrix() {
  local idx_dir="$1"
  local outfile="$2"
  local length_file="${3:-}"
  local -a files=()
  shopt -s nullglob
  files=("${idx_dir}"/*.idxstats.tsv)
  shopt -u nullglob
  if (( ${#files[@]} == 0 )); then
    echo "ERROR: no idxstats files in $idx_dir" >&2
    exit 1
  fi
  awk -v length_file="$length_file" '
    FNR == 1 {
      nfiles++
      fname = FILENAME
      sub(/^.*\//, "", fname)
      sub(/\.idxstats\.tsv$/, "", fname)
      samples[nfiles] = fname
    }
    $1 != "*" {
      key = $1
      if (!(key in seen)) {
        nkeys++
        order[nkeys] = key
        seen[key] = 1
        lengths[key] = $2
      }
      counts[key, nfiles] = $3
    }
    END {
      printf "transcript"
      for (i = 1; i <= nfiles; i++) printf "\t%s", samples[i]
      printf "\n"
      for (k = 1; k <= nkeys; k++) {
        key = order[k]
        printf "%s", key
        for (i = 1; i <= nfiles; i++) printf "\t%d", counts[key, i] + 0
        printf "\n"
      }
      if (length_file != "") {
        print "transcript\tlength" > length_file
        for (k = 1; k <= nkeys; k++) {
          key = order[k]
          printf "%s\t%d\n", key, lengths[key] + 0 > length_file
        }
        close(length_file)
      }
    }
  ' "${files[@]}" > "$outfile"
}

write_size_factors() {
  local y_table="$1"
  local outfile="$2"
  local min_y="$3"
  awk -v min_y="$min_y" 'BEGIN { OFS="\t" }
    {
      sample[++n] = $1
      assay[n] = $2
      y[n] = $3 + 0
      if (y[n] > 0) { slog += log(y[n]); npos++ }
      if (y[n] < min_y) warn[++nw] = $1 " (" $2 ")\t" y[n]
    }
    END {
      if (n == 0) {
        print "ERROR: no yeast counts" > "/dev/stderr"
        exit 1
      }
      if (npos == 0) {
        print "ERROR: all yeast primary-mapped counts are 0" > "/dev/stderr"
        exit 1
      }
      gm = exp(slog / npos)
      print "sample", "assay", "Y", "size_factor"
      for (i = 1; i <= n; i++) {
        sf = (y[i] > 0 ? sprintf("%.6f", y[i] / gm) : "NA")
        print sample[i], assay[i], y[i], sf
      }
      if (nw) {
        print "WARNING: libraries below MIN_Y_READS=" min_y ":" > "/dev/stderr"
        for (i = 1; i <= nw; i++) print "  " warn[i] > "/dev/stderr"
      }
    }
  ' "$y_table" > "$outfile"
}

RIBO_QC_DIR="${OUTPUT_DIR}/01.clean_fastq/ribo"
RIBO_SMRNA_DIR="${OUTPUT_DIR}/02.no_smrna/ribo"
RIBO_MAP_DIR="${OUTPUT_DIR}/03.mapped/ribo"
RIBO_FILTER_DIR="${OUTPUT_DIR}/04.filtered_bam/ribo"
RIBO_COUNT_DIR="${OUTPUT_DIR}/05.count/ribo"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$RIBO_QC_DIR" "$RIBO_SMRNA_DIR" "$RIBO_MAP_DIR" \
  "$RIBO_FILTER_DIR" "$RIBO_COUNT_DIR" "$LOG_DIR"

if [[ "$WITH_RNA" == "1" ]]; then
  RNA_QC_DIR="${OUTPUT_DIR}/01.clean_fastq/rna"
  RNA_MAP_DIR="${OUTPUT_DIR}/03.mapped/rna"
  RNA_FILTER_DIR="${OUTPUT_DIR}/04.filtered_bam/rna"
  RNA_COUNT_DIR="${OUTPUT_DIR}/05.count/rna"
  mkdir -p "$RNA_QC_DIR" "$RNA_MAP_DIR" "$RNA_FILTER_DIR" "$RNA_COUNT_DIR"
fi
if [[ "$SPIKEIN" == "1" ]]; then
  YEAST_MAP_RIBO_DIR="${OUTPUT_DIR}/03.mapped_yeast/ribo"
  YEAST_FILTER_RIBO_DIR="${OUTPUT_DIR}/04.filtered_bam_yeast/ribo"
  SPIKEIN_DIR="${OUTPUT_DIR}/06.spikein"
  mkdir -p "$YEAST_MAP_RIBO_DIR" "$YEAST_FILTER_RIBO_DIR" "$SPIKEIN_DIR"
  if [[ "$WITH_RNA" == "1" ]]; then
    YEAST_MAP_RNA_DIR="${OUTPUT_DIR}/03.mapped_yeast/rna"
    YEAST_FILTER_RNA_DIR="${OUTPUT_DIR}/04.filtered_bam_yeast/rna"
    mkdir -p "$YEAST_MAP_RNA_DIR" "$YEAST_FILTER_RNA_DIR"
  fi
fi

mapfile -t ribo_r1_files < <(collect_r1_files "$RIBO_RAW_DIR")
echo "Found ${#ribo_r1_files[@]} Ribo-seq paired-end sample(s)."

for r1 in "${ribo_r1_files[@]}"; do
  sample="$(sample_name_from_r1 "$r1")"
  if ! r2="$(find_r2 "$r1")"; then
    echo "ERROR: R2 file missing for Ribo-seq sample $sample" >&2
    exit 1
  fi

  clean_fq="${RIBO_QC_DIR}/${sample}.fastq.gz"
  nosmrna_fq="${RIBO_SMRNA_DIR}/${sample}.fastq.gz"
  mapped_bam="${RIBO_MAP_DIR}/${sample}.bam"
  filtered_bam="${RIBO_FILTER_DIR}/${sample}.bam"
  idxstats="${RIBO_COUNT_DIR}/${sample}.idxstats.tsv"

  echo "[$sample] 1/5 fastp merge + seqkit length filter"
  fastp \
    --in1 "$r1" \
    --in2 "$r2" \
    --out1 /dev/null \
    --out2 /dev/null \
    --n_base_limit 5 \
    --cut_right \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 25 \
    --trim_poly_g \
    --poly_x_min_len 15 \
    --length_required "$MIN_LENGTH" \
    --thread "$THREADS" \
    --html "${RIBO_QC_DIR}/${sample}.fastp.html" \
    --json "${RIBO_QC_DIR}/${sample}.fastp.json" \
    --overlap_len_require 15 \
    --merge \
    --merged_out /dev/stdout \
    2> "${LOG_DIR}/${sample}.fastp.log" |
  seqkit seq \
    -t dna \
    -m "$MIN_LENGTH" \
    -M "$MAX_LENGTH" \
    -g \
    -j "$THREADS" \
    -o "$clean_fq" -

  echo "[$sample] 2/5 Bowtie rRNA/tRNA/snRNA depletion"
  bowtie \
    -q \
    -v "$BOWTIE_MISMATCHES" \
    -p "$THREADS" \
    -m "$BOWTIE_MAX_ALIGNMENTS" \
    --best --strata \
    -x "$SMRNA_INDEX" \
    --un "${RIBO_SMRNA_DIR}/${sample}.unmapped.fastq" \
    -S "$clean_fq" \
    2> "${LOG_DIR}/${sample}.bowtie_smrna.log" |
  samtools view -@ "$THREADS" -b -F 4 -o "${RIBO_SMRNA_DIR}/${sample}.smrna.bam" -
  gzip -f "${RIBO_SMRNA_DIR}/${sample}.unmapped.fastq"
  mv -f "${RIBO_SMRNA_DIR}/${sample}.unmapped.fastq.gz" "$nosmrna_fq"

  echo "[$sample] 3/5 HISAT2 transcriptome mapping"
  hisat2_cmd=(
    hisat2
    -p "$THREADS"
    --rg-id "$sample"
    --rg "SM:${sample}"
    --rg "LB:${sample}"
    --rg "PL:ILLUMINA"
    -x "$TX_INDEX"
    -U "$nosmrna_fq"
  )
  if [[ -n "$HISAT2_STRANDNESS" ]]; then
    hisat2_cmd+=(--rna-strandness "$HISAT2_STRANDNESS")
  fi
  if [[ "$SPIKEIN" == "1" ]]; then
    hisat2_cmd+=(--un-gz "${RIBO_MAP_DIR}/${sample}_unmapped.fastq.gz")
  fi
  "${hisat2_cmd[@]}" \
    2> "${LOG_DIR}/${sample}.hisat2.log" |
  samtools sort -@ "$THREADS" -O BAM -o "$mapped_bam" -
  samtools index -@ "$THREADS" "$mapped_bam"

  echo "[$sample] 4/5 samtools filtering"
  filter_sorted_bam "$mapped_bam" "$filtered_bam" "$MIN_LENGTH" "$MAX_LENGTH"

  echo "[$sample] 5/5 samtools idxstats"
  samtools idxstats "$filtered_bam" > "$idxstats"
done

merge_idxstats_matrix "$RIBO_COUNT_DIR" \
  "${RIBO_COUNT_DIR}/transcript_count_matrix.txt" \
  "${RIBO_COUNT_DIR}/transcript_lengths.tsv"
echo "Ribo-seq transcript counts: ${RIBO_COUNT_DIR}/transcript_count_matrix.txt"

if [[ "$WITH_RNA" == "1" ]]; then
  mapfile -t rna_r1_files < <(collect_r1_files "$RNA_RAW_DIR")
  echo "Found ${#rna_r1_files[@]} RNA-seq paired-end sample(s)."

  for r1 in "${rna_r1_files[@]}"; do
    sample="$(sample_name_from_r1 "$r1")"
    if ! r2="$(find_r2 "$r1")"; then
      echo "ERROR: R2 file missing for RNA-seq sample $sample" >&2
      exit 1
    fi

    clean_r1="${RNA_QC_DIR}/${sample}.R1.fastq.gz"
    clean_r2="${RNA_QC_DIR}/${sample}.R2.fastq.gz"
    mapped_bam="${RNA_MAP_DIR}/${sample}.bam"
    filtered_bam="${RNA_FILTER_DIR}/${sample}.bam"
    idxstats="${RNA_COUNT_DIR}/${sample}.idxstats.tsv"

    echo "[$sample] 1/4 fastp paired-end QC"
    fastp \
      --in1 "$r1" \
      --in2 "$r2" \
      --out1 "$clean_r1" \
      --out2 "$clean_r2" \
      --n_base_limit 5 \
      --cut_right \
      --cut_tail \
      --cut_window_size 4 \
      --cut_mean_quality 25 \
      --trim_poly_g \
      --poly_x_min_len 15 \
      --length_required "$RNA_MIN_LENGTH" \
      --length_limit "$RNA_MAX_LENGTH" \
      --thread "$THREADS" \
      --html "${RNA_QC_DIR}/${sample}.fastp.html" \
      --json "${RNA_QC_DIR}/${sample}.fastp.json" \
      2> "${LOG_DIR}/${sample}.fastp.log"

    echo "[$sample] 2/4 HISAT2 transcriptome mapping"
    hisat2_cmd=(
      hisat2
      -p "$THREADS"
      --rg-id "$sample"
      --rg "SM:${sample}"
      --rg "LB:${sample}"
      --rg "PL:ILLUMINA"
      -x "$TX_INDEX"
      -1 "$clean_r1"
      -2 "$clean_r2"
    )
    if [[ -n "$HISAT2_RNA_STRANDNESS" ]]; then
      hisat2_cmd+=(--rna-strandness "$HISAT2_RNA_STRANDNESS")
    fi
    "${hisat2_cmd[@]}" \
      2> "${LOG_DIR}/${sample}.hisat2.log" |
    samtools sort -@ "$THREADS" -O BAM -o "$mapped_bam" -
    samtools index -@ "$THREADS" "$mapped_bam"

    echo "[$sample] 3/4 samtools filtering"
    filter_sorted_bam "$mapped_bam" "$filtered_bam"

    echo "[$sample] 4/4 samtools idxstats"
    samtools idxstats "$filtered_bam" > "$idxstats"
  done

  merge_idxstats_matrix "$RNA_COUNT_DIR" \
    "${RNA_COUNT_DIR}/transcript_count_matrix.txt" \
    "${RNA_COUNT_DIR}/transcript_lengths.tsv"
  echo "RNA-seq transcript counts: ${RNA_COUNT_DIR}/transcript_count_matrix.txt"
fi

if [[ "$SPIKEIN" == "1" ]]; then
  y_table="${SPIKEIN_DIR}/yeast_primary_counts.tsv"
  : > "$y_table"

  shopt -s nullglob
  ribo_unmapped=("${RIBO_MAP_DIR}"/*_unmapped.fastq.gz)
  shopt -u nullglob
  if (( ${#ribo_unmapped[@]} == 0 )); then
    echo "ERROR: no host-unmapped Ribo FASTQ files in $RIBO_MAP_DIR" >&2
    exit 1
  fi

  echo "Mapping host-unmapped Ribo reads to yeast transcriptome"
  for fq in "${ribo_unmapped[@]}"; do
    sample="$(basename -- "$fq")"
    sample="${sample%_unmapped.fastq.gz}"
    yeast_bam="${YEAST_MAP_RIBO_DIR}/${sample}.bam"
    yeast_filtered="${YEAST_FILTER_RIBO_DIR}/${sample}.bam"

    hisat2_cmd=(
      hisat2
      -p "$THREADS"
      --rg-id "$sample"
      --rg "SM:${sample}"
      --rg "LB:${sample}"
      --rg "PL:ILLUMINA"
      -x "$YEAST_TX_INDEX"
      -U "$fq"
    )
    if [[ -n "$YEAST_HISAT2_STRANDNESS" ]]; then
      hisat2_cmd+=(--rna-strandness "$YEAST_HISAT2_STRANDNESS")
    fi
    "${hisat2_cmd[@]}" \
      2> "${LOG_DIR}/${sample}.yeast_hisat2.log" |
    samtools sort -@ "$THREADS" -O BAM -o "$yeast_bam" -
    samtools index -@ "$THREADS" "$yeast_bam"
    filter_sorted_bam "$yeast_bam" "$yeast_filtered"

    y="$(samtools view -c -q "$MIN_MAPQ" -F 2308 "$yeast_filtered")"
    printf '%s\tribo\t%s\n' "$sample" "$y" >> "$y_table"
  done

  if [[ "$WITH_RNA" == "1" ]]; then
    echo "Mapping clean RNA reads to yeast transcriptome"
    shopt -s nullglob
    rna_clean_r1=("${RNA_QC_DIR}"/*.R1.fastq.gz)
    shopt -u nullglob
    for r1 in "${rna_clean_r1[@]}"; do
      sample="$(sample_name_from_r1 "$r1")"
      r2="${RNA_QC_DIR}/${sample}.R2.fastq.gz"
      yeast_bam="${YEAST_MAP_RNA_DIR}/${sample}.bam"
      yeast_filtered="${YEAST_FILTER_RNA_DIR}/${sample}.bam"

      hisat2_cmd=(
        hisat2
        -p "$THREADS"
        --rg-id "$sample"
        --rg "SM:${sample}"
        --rg "LB:${sample}"
        --rg "PL:ILLUMINA"
        -x "$YEAST_TX_INDEX"
        -1 "$r1"
        -2 "$r2"
      )
      if [[ -n "$YEAST_HISAT2_RNA_STRANDNESS" ]]; then
        hisat2_cmd+=(--rna-strandness "$YEAST_HISAT2_RNA_STRANDNESS")
      fi
      "${hisat2_cmd[@]}" \
        2> "${LOG_DIR}/${sample}.yeast_hisat2.log" |
      samtools sort -@ "$THREADS" -O BAM -o "$yeast_bam" -
      samtools index -@ "$THREADS" "$yeast_bam"
      filter_sorted_bam "$yeast_bam" "$yeast_filtered"

      y="$(samtools view -c -q "$MIN_MAPQ" -F 2308 "$yeast_filtered")"
      printf '%s\trna\t%s\n' "$sample" "$y" >> "$y_table"
    done
  fi

  write_size_factors "$y_table" "${SPIKEIN_DIR}/spikein_size_factors.tsv" "$MIN_Y_READS"
  echo "Spike-in size factors: ${SPIKEIN_DIR}/spikein_size_factors.tsv"
fi

if [[ "$WITH_RNA" == "1" ]]; then
  if [[ -z "$SAMPLE_SHEET" ]]; then
    echo "WARNING: SAMPLE_SHEET is empty; skipping manual TE. Copy ribo-seq.samples.example.tsv." >&2
  elif [[ ! -f "$SAMPLE_SHEET" ]]; then
    echo "ERROR: SAMPLE_SHEET not found: $SAMPLE_SHEET" >&2
    exit 1
  elif ! command -v "$RSCRIPT" >/dev/null 2>&1; then
    echo "ERROR: Rscript is not available in PATH (needed for manual TE): $RSCRIPT" >&2
    exit 1
  else
    TE_DIR="${OUTPUT_DIR}/07.TE"
    te_cmd=(
      "$RSCRIPT" "${SCRIPT_DIR}/calc_te.R"
      --mode ribo
      --ribo-counts "${RIBO_COUNT_DIR}/transcript_count_matrix.txt"
      --rna-counts "${RNA_COUNT_DIR}/transcript_count_matrix.txt"
      --lengths "${RIBO_COUNT_DIR}/transcript_lengths.tsv"
      --sample-sheet "$SAMPLE_SHEET"
      --output-dir "$TE_DIR"
      --pseudocount "$TE_PSEUDOCOUNT"
      --min-tpm "$TE_MIN_TPM"
    )
    if [[ "$SPIKEIN" == "1" && -f "${SPIKEIN_DIR}/spikein_size_factors.tsv" ]]; then
      te_cmd+=(--spikein-table "${SPIKEIN_DIR}/spikein_size_factors.tsv")
    fi
    echo "Calculating manual TE"
    "${te_cmd[@]}"
    echo "Manual TE tables: $TE_DIR"
  fi
fi

echo "Pipeline completed. Count tables: ${OUTPUT_DIR}/05.count"
