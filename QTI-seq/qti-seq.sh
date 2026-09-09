#!/usr/bin/env bash
# QTI-seq upstream workflow: raw paired-end FASTQ to genome CDS count matrices.
# QC (fastp merge) -> rRNA/tRNA depletion (Bowtie) -> STAR unique genome mapping
# -> MAPQ/flag/length filter (samtools) -> CDS counts (featureCounts).
# Optional matched RNA-seq/Input libraries are kept paired-end.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/qti-seq.config.sh}"

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: configuration file not found: $CONFIG" >&2
  echo "Copy qti-seq.config.example.sh to qti-seq.config.sh and edit the paths first." >&2
  exit 1
fi

# shellcheck source=/dev/null
source "$CONFIG"

WITH_RNA="${WITH_RNA:-0}"

required_variables=(
  QTI_RAW_DIR OUTPUT_DIR SMRNA_INDEX STAR_INDEX GTF THREADS
  STAR_SJDB_OVERHANG STAR_MISMATCHES STAR_BAM_SORT_RAM
  MIN_LENGTH MAX_LENGTH BOWTIE_MISMATCHES BOWTIE_MAX_ALIGNMENTS
  MIN_MAPQ FILTER_FLAGS FC_FEATURE FC_ATTRIBUTE FC_STRANDNESS
)
if [[ "$WITH_RNA" == "1" ]]; then
  required_variables+=(RNA_RAW_DIR RNA_MIN_LENGTH RNA_MAX_LENGTH)
fi
for variable in "${required_variables[@]}"; do
  if [[ -z "${!variable:-}" ]]; then
    echo "ERROR: required configuration variable is empty: $variable" >&2
    exit 1
  fi
done

for program in fastp seqkit bowtie STAR samtools gzip featureCounts; do
  if ! command -v "$program" >/dev/null 2>&1; then
    echo "ERROR: required program is not available in PATH: $program" >&2
    exit 1
  fi
done

if [[ ! -d "$QTI_RAW_DIR" ]]; then
  echo "ERROR: QTI-seq FASTQ directory not found: $QTI_RAW_DIR" >&2
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
if [[ ! -f "${STAR_INDEX}/Genome" ]]; then
  echo "ERROR: STAR genome index not found: ${STAR_INDEX}/Genome" >&2
  exit 1
fi
if [[ ! -f "$GTF" ]]; then
  echo "ERROR: GTF annotation not found: $GTF" >&2
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
  local dir sample candidate
  dir="$(dirname -- "$r1")"
  sample="$(sample_name_from_r1 "$r1")"
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

run_star() {
  local prefix="$1"
  local out_dir="$2"
  shift 2
  local -a reads=("$@")
  local star_prefix="${out_dir}/${prefix}_"
  STAR \
    --runThreadN "$THREADS" \
    --runMode alignReads \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "${reads[@]}" \
    --readFilesCommand zcat \
    --outFilterMismatchNmax "$STAR_MISMATCHES" \
    --outFilterMultimapNmax 1 \
    --outSAMattributes All \
    --outSAMattrRGline "ID:${prefix} SM:${prefix} LB:${prefix} PL:ILLUMINA" \
    --alignEndsType EndToEnd \
    --sjdbGTFfile "$GTF" \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --alignIntronMax 2000000 \
    --sjdbOverhang "$STAR_SJDB_OVERHANG" \
    --outSJfilterReads Unique \
    --outReadsUnmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM "$STAR_BAM_SORT_RAM" \
    --outMultimapperOrder Random \
    --outSAMmultNmax 1 \
    --outFileNamePrefix "$star_prefix" \
    > "${LOG_DIR}/${prefix}.star.log" 2>&1
  samtools index -@ "$THREADS" "${star_prefix}Aligned.sortedByCoord.out.bam"
  rm -rf "${out_dir}/${prefix}__STARtmp"
}

QTI_QC_DIR="${OUTPUT_DIR}/01.clean_fastq/qti"
QTI_SMRNA_DIR="${OUTPUT_DIR}/02.no_smrna/qti"
QTI_MAP_DIR="${OUTPUT_DIR}/03.mapped/qti"
QTI_FILTER_DIR="${OUTPUT_DIR}/04.filtered_bam/qti"
COUNT_DIR="${OUTPUT_DIR}/05.count"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$QTI_QC_DIR" "$QTI_SMRNA_DIR" "$QTI_MAP_DIR" \
  "$QTI_FILTER_DIR" "$COUNT_DIR" "$LOG_DIR"

if [[ "$WITH_RNA" == "1" ]]; then
  RNA_QC_DIR="${OUTPUT_DIR}/01.clean_fastq/rna"
  RNA_MAP_DIR="${OUTPUT_DIR}/03.mapped/rna"
  RNA_FILTER_DIR="${OUTPUT_DIR}/04.filtered_bam/rna"
  mkdir -p "$RNA_QC_DIR" "$RNA_MAP_DIR" "$RNA_FILTER_DIR"
fi

mapfile -t qti_r1_files < <(collect_r1_files "$QTI_RAW_DIR")
echo "Found ${#qti_r1_files[@]} QTI-seq paired-end sample(s)."

for r1 in "${qti_r1_files[@]}"; do
  sample="$(sample_name_from_r1 "$r1")"
  if ! r2="$(find_r2 "$r1")"; then
    echo "ERROR: R2 file missing for QTI-seq sample $sample" >&2
    exit 1
  fi

  clean_fq="${QTI_QC_DIR}/${sample}.fastq.gz"
  nosmrna_fq="${QTI_SMRNA_DIR}/${sample}.fastq.gz"
  mapped_bam="${QTI_MAP_DIR}/${sample}_Aligned.sortedByCoord.out.bam"
  filtered_bam="${QTI_FILTER_DIR}/${sample}.bam"

  echo "[$sample] 1/4 fastp merge + seqkit length filter"
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
    --html "${QTI_QC_DIR}/${sample}.fastp.html" \
    --json "${QTI_QC_DIR}/${sample}.fastp.json" \
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

  echo "[$sample] 2/4 Bowtie rRNA/tRNA/snRNA depletion"
  bowtie \
    -q \
    -v "$BOWTIE_MISMATCHES" \
    -p "$THREADS" \
    -m "$BOWTIE_MAX_ALIGNMENTS" \
    --best --strata \
    -x "$SMRNA_INDEX" \
    --un "${QTI_SMRNA_DIR}/${sample}.unmapped.fastq" \
    -S "$clean_fq" \
    2> "${LOG_DIR}/${sample}.bowtie_smrna.log" |
  samtools view -@ "$THREADS" -b -F 4 -o "${QTI_SMRNA_DIR}/${sample}.smrna.bam" -
  gzip -f "${QTI_SMRNA_DIR}/${sample}.unmapped.fastq"
  mv -f "${QTI_SMRNA_DIR}/${sample}.unmapped.fastq.gz" "$nosmrna_fq"

  echo "[$sample] 3/4 STAR unique genome mapping"
  run_star "$sample" "$QTI_MAP_DIR" "$nosmrna_fq"

  echo "[$sample] 4/4 samtools filtering"
  filter_sorted_bam "$mapped_bam" "$filtered_bam" "$MIN_LENGTH" "$MAX_LENGTH"
done

shopt -s nullglob
qti_bams=("${QTI_FILTER_DIR}"/*.bam)
shopt -u nullglob
if (( ${#qti_bams[@]} == 0 )); then
  echo "ERROR: no filtered QTI BAM files in $QTI_FILTER_DIR" >&2
  exit 1
fi

echo "Counting QTI-seq CDS reads with featureCounts"
featureCounts \
  -T "$THREADS" \
  -s "$FC_STRANDNESS" \
  -t "$FC_FEATURE" \
  -g "$FC_ATTRIBUTE" \
  -a "$GTF" \
  -o "${COUNT_DIR}/qti_cds_counts.txt" \
  "${qti_bams[@]}" \
  2> "${LOG_DIR}/featureCounts.qti.log"

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
    mapped_bam="${RNA_MAP_DIR}/${sample}_Aligned.sortedByCoord.out.bam"
    filtered_bam="${RNA_FILTER_DIR}/${sample}.bam"

    echo "[$sample] 1/3 fastp paired-end QC"
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

    echo "[$sample] 2/3 STAR unique genome mapping"
    run_star "$sample" "$RNA_MAP_DIR" "$clean_r1" "$clean_r2"

    echo "[$sample] 3/3 samtools filtering"
    filter_sorted_bam "$mapped_bam" "$filtered_bam"
  done

  shopt -s nullglob
  rna_bams=("${RNA_FILTER_DIR}"/*.bam)
  shopt -u nullglob
  if (( ${#rna_bams[@]} == 0 )); then
    echo "ERROR: no filtered RNA BAM files in $RNA_FILTER_DIR" >&2
    exit 1
  fi

  echo "Counting RNA-seq CDS fragments with featureCounts"
  featureCounts \
    -T "$THREADS" \
    -s "$FC_STRANDNESS" \
    -p --countReadPairs \
    -t "$FC_FEATURE" \
    -g "$FC_ATTRIBUTE" \
    -a "$GTF" \
    -o "${COUNT_DIR}/rna_cds_counts.txt" \
    "${rna_bams[@]}" \
    2> "${LOG_DIR}/featureCounts.rna.log"
fi

echo "Pipeline completed. Count tables: $COUNT_DIR"
