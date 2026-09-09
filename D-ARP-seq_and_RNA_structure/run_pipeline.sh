#!/usr/bin/env bash
# DHU-seq upstream workflow: raw paired-end FASTQ to merged PMSM tables.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/config.sh}"

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: configuration file not found: $CONFIG" >&2
  echo "Copy config.example.sh to config.sh and edit the paths first." >&2
  exit 1
fi

# shellcheck source=/dev/null
source "$CONFIG"

required_variables=(
  RAW_DIR OUTPUT_DIR BOWTIE_INDEX REFERENCE_FASTA SAMPLE_SHEET THREADS JAVA_MEMORY
  FIRST_PASS_MIN_LENGTH MIN_LENGTH MAX_LENGTH UMI_LOCATION UMI_LENGTH UMI_PREFIX
  BOWTIE_MISMATCHES BOWTIE_MAX_ALIGNMENTS MIN_MAPQ FILTER_FLAGS
  UMI_EDIT_DISTANCE UMI_ALGORITHM UMI_SEPARATOR
)
for variable in "${required_variables[@]}"; do
  if [[ -z "${!variable:-}" ]]; then
    echo "ERROR: required configuration variable is empty: $variable" >&2
    exit 1
  fi
done

for program in fastp seqkit bowtie bowtie-build samtools bedtools umicollapse python Rscript awk sha256sum; do
  if ! command -v "$program" >/dev/null 2>&1; then
    echo "ERROR: required program is not available in PATH: $program" >&2
    exit 1
  fi
done

if [[ ! -d "$RAW_DIR" ]]; then
  echo "ERROR: FASTQ directory not found: $RAW_DIR" >&2
  exit 1
fi
if [[ ! -f "$REFERENCE_FASTA" ]]; then
  echo "ERROR: reference FASTA not found: $REFERENCE_FASTA" >&2
  exit 1
fi
if [[ ! -f "$SAMPLE_SHEET" ]]; then
  echo "ERROR: sample sheet not found: $SAMPLE_SHEET" >&2
  exit 1
fi
if [[ -n "${EXPECTED_REFERENCE_SHA256:-}" ]]; then
  actual_reference_sha256="$(sha256sum "$REFERENCE_FASTA" | awk '{print $1}')"
  if [[ "$actual_reference_sha256" != "$EXPECTED_REFERENCE_SHA256" ]]; then
    echo "ERROR: reference FASTA checksum does not match." >&2
    echo "  expected: $EXPECTED_REFERENCE_SHA256" >&2
    echo "  observed: $actual_reference_sha256" >&2
    exit 1
  fi
fi
if [[ ! -f "${BOWTIE_INDEX}.1.ebwt" && ! -f "${BOWTIE_INDEX}.1.ebwtl" ]]; then
  echo "Bowtie 1 index is absent; building it from $REFERENCE_FASTA"
  mkdir -p "$(dirname -- "$BOWTIE_INDEX")"
  bowtie-build "$REFERENCE_FASTA" "$BOWTIE_INDEX"
fi

if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
  echo "Indexing reference FASTA with samtools faidx"
  samtools faidx "$REFERENCE_FASTA"
fi
REFERENCE_FAI="${REFERENCE_FASTA}.fai"

RMUMI_DIR="${OUTPUT_DIR}/2.rmumi_fq"
QC_DIR="${OUTPUT_DIR}/3.clean_fq"
MAP_DIR="${OUTPUT_DIR}/4.mapped"
FILTER_DIR="${OUTPUT_DIR}/5.filtered_bam"
DEDUP_DIR="${OUTPUT_DIR}/6.bam_rmdup"
STOP_DIR="${OUTPUT_DIR}/8.stop"
MUTATION_DIR="${OUTPUT_DIR}/9.mutation"
PMSM_DIR="${OUTPUT_DIR}/10.merge"
ANNO_DIR="${OUTPUT_DIR}/11.anno_sprinzl"
BACKGROUND_DIR="${OUTPUT_DIR}/13_calu_bg_tRNA"
MERGE_REP_DIR="${OUTPUT_DIR}/14_merge_DHU_site_rep"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$RMUMI_DIR" "$QC_DIR" "$MAP_DIR" "$FILTER_DIR" "$DEDUP_DIR" \
  "$STOP_DIR" "$MUTATION_DIR" "$PMSM_DIR" "$ANNO_DIR" \
  "$BACKGROUND_DIR" "$MERGE_REP_DIR" "$LOG_DIR"

mapfile -t sample_ids < <(
  awk -F, 'NR > 1 && $1 != "" { sub(/\r$/, "", $1); print $1 }' "$SAMPLE_SHEET" |
    sort -u
)
if (( ${#sample_ids[@]} == 0 )); then
  echo "ERROR: sample sheet contains no sample_id values: $SAMPLE_SHEET" >&2
  exit 1
fi

r1_files=()
for sample in "${sample_ids[@]}"; do
  mapfile -t matches < <(
    find "$RAW_DIR" \( -type f -o -type l \) -name "${sample}.R1.fastq.gz" -print
  )
  if (( ${#matches[@]} != 1 )); then
    echo "ERROR: expected exactly one R1 for $sample under $RAW_DIR; found ${#matches[@]}" >&2
    printf '  %s\n' "${matches[@]}" >&2
    exit 1
  fi
  r1_files+=("${matches[0]}")
done

echo "Found ${#r1_files[@]} paired-end sample(s)."
for r1 in "${r1_files[@]}"; do
  filename="$(basename -- "$r1")"
  sample="${filename%.R1.fastq.gz}"
  r2="$(dirname -- "$r1")/${sample}.R2.fastq.gz"
  if [[ ! -f "$r2" ]]; then
    echo "ERROR: R2 file missing for sample $sample: $r2" >&2
    exit 1
  fi

  rmumi_fq="${RMUMI_DIR}/${sample}-merged-filtered-revcomp.R1.fastq.gz"
  clean_fq="${QC_DIR}/${sample}-merged-filtered-revcomp-filtered.R1.fastq.gz"
  pipeline_prefix="${sample}-merged-filtered-revcomp-filtered.R1"
  mapped_bam="${MAP_DIR}/${pipeline_prefix}.bam"
  filtered_bam="${FILTER_DIR}/${pipeline_prefix}.bam"
  dedup_bam="${DEDUP_DIR}/${pipeline_prefix}-dedup.bam"
  fastp_unmerged_r1="${RMUMI_DIR}/${sample}.R1.fastq.gz"
  fastp_unmerged_r2="${RMUMI_DIR}/${sample}.R2.fastq.gz"

  echo "[$sample] 1/5 first QC: UMI + paired-end merge + reverse complement"
  fastp \
    --in1 "$r1" \
    --in2 "$r2" \
    --out1 "$fastp_unmerged_r1" \
    --out2 "$fastp_unmerged_r2" \
    --umi \
    --umi_loc "$UMI_LOCATION" \
    --umi_len "$UMI_LENGTH" \
    --umi_prefix "$UMI_PREFIX" \
    --n_base_limit 5 \
    --cut_right \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 25 \
    --trim_poly_g \
    --poly_x_min_len 15 \
    --length_required "$FIRST_PASS_MIN_LENGTH" \
    --thread "$THREADS" \
    --html "${RMUMI_DIR}/${sample}_clean.html" \
    --json "${RMUMI_DIR}/${sample}.fastp.json" \
    --overlap_len_require 15 \
    --merge \
    --merged_out /dev/stdout \
    2> "${LOG_DIR}/${sample}.fastp_pass1.log" |
  seqkit seq \
    -t dna \
    -m "$FIRST_PASS_MIN_LENGTH" \
    -M "$MAX_LENGTH" \
    -g -p -r \
    -j "$THREADS" \
    -o "$rmumi_fq" -
  rm -f -- "$fastp_unmerged_r1" "$fastp_unmerged_r2"

  echo "[$sample] 2/5 second QC: quality and 20-200 nt filter"
  fastp \
    --in1 "$rmumi_fq" \
    --out1 "$clean_fq" \
    --n_base_limit 5 \
    --cut_right \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 25 \
    --trim_poly_g \
    --poly_x_min_len 15 \
    --length_required "$MIN_LENGTH" \
    --length_limit "$MAX_LENGTH" \
    --thread "$THREADS" \
    --html "${QC_DIR}/${sample}-merged-filtered-revcomp_clean.html" \
    --json "${QC_DIR}/${sample}-merged-filtered-revcomp.fastp.json" \
    2> "${LOG_DIR}/${sample}.fastp_pass2.log"

  echo "[$sample] 3/5 Bowtie mapping"
  bowtie \
    -p "$THREADS" \
    -q -S \
    --sam-RG "ID:${pipeline_prefix}" \
    --sam-RG "SM:${pipeline_prefix}" \
    --sam-RG "LB:${pipeline_prefix}" \
    --sam-RG "PL:ILLUMINA" \
    -v "$BOWTIE_MISMATCHES" \
    -m "$BOWTIE_MAX_ALIGNMENTS" \
    --best --strata --norc \
    "$BOWTIE_INDEX" "$clean_fq" \
    2> "${LOG_DIR}/${sample}.bowtie.log" |
  samtools sort \
    -@ "$THREADS" \
    -O BAM \
    -o "$mapped_bam" -
  samtools index -@ "$THREADS" "$mapped_bam"

  echo "[$sample] 4/5 samtools filtering"
  samtools view \
    -@ "$THREADS" \
    -h \
    -q "$MIN_MAPQ" \
    -F "$FILTER_FLAGS" \
    "$mapped_bam" |
  awk 'BEGIN { OFS="\t" } /^@/ || $6 !~ /[ID]/' |
  samtools view \
    -@ "$THREADS" \
    -b \
    -o "$filtered_bam" -
  samtools index -@ "$THREADS" "$filtered_bam"

  echo "[$sample] 5/5 UMICollapse deduplication"
  env _JAVA_OPTIONS="-Xmx${JAVA_MEMORY} -Xms512M -Xss2M -XX:+UseG1GC -XX:NewRatio=1 -XX:+UseStringDeduplication" \
    umicollapse bam \
      -i "$filtered_bam" \
      -o "$dedup_bam" \
      -k "$UMI_EDIT_DISTANCE" \
      --algo "$UMI_ALGORITHM" \
      --umi-sep "$UMI_SEPARATOR" \
      -u "$UMI_LENGTH"
  samtools index -@ "$THREADS" "$dedup_bam"
done

echo "Generating RT-stop signal with the original call_stop.py"
python "${SCRIPT_DIR}/pmsm/call_stop.py" \
  --input-dir "$DEDUP_DIR" \
  --output-dir "$STOP_DIR" \
  --genome-file "$REFERENCE_FAI" \
  --end-type 5 \
  --plus-shift 0 \
  --minus-shift 0 \
  --workers "$THREADS"

echo "Generating mutation/PMSM components with the original mutation_FR_v4.py"
python "${SCRIPT_DIR}/pmsm/mutation_FR_v4.py" \
  --input-dir "$DEDUP_DIR" \
  --output-dir "$MUTATION_DIR" \
  --log-dir "${MUTATION_DIR}/logs" \
  --reference "$REFERENCE_FASTA" \
  --parallel-jobs "$THREADS" \
  --base-quality 0 \
  --mapping-quality 0

echo "Merging stop and mutation signals with the original merge_stop_mutation_v2.py"
python "${SCRIPT_DIR}/pmsm/merge_stop_mutation_v2.py" \
  --stop-dir "$STOP_DIR" \
  --mutation-dir "$MUTATION_DIR" \
  --output-dir "$PMSM_DIR" \
  --fai-file "$REFERENCE_FAI"

echo "Adding Sprinzl positions with the original annotation code"
python "${SCRIPT_DIR}/postprocess/annotate_tRNA_position_sprinzl.py" \
  --input "$PMSM_DIR" \
  --output "$ANNO_DIR" \
  --map-file "${SCRIPT_DIR}/postprocess/hg38_trna_raw_to_sprinzl.tsv" \
  --map-id-col reference_id \
  --pattern "*.txt" \
  --workers "$THREADS"

echo "Applying background correction and merging biological replicates"
Rscript "${SCRIPT_DIR}/postprocess/run_postprocess.R" \
  "$SAMPLE_SHEET" \
  "$ANNO_DIR" \
  "$BACKGROUND_DIR" \
  "$MERGE_REP_DIR"

echo "Pipeline completed. Replicate-merged tables: $MERGE_REP_DIR"
