#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: Rscript run_postprocess.R SAMPLE_SHEET ANNOTATED_DIR ",
    "BACKGROUND_OUTPUT_DIR MERGED_OUTPUT_DIR"
  )
}

sample_sheet_path <- normalizePath(args[[1]], mustWork = TRUE)
annotated_dir <- normalizePath(args[[2]], mustWork = TRUE)
background_output_dir <- args[[3]]
merged_output_dir <- args[[4]]

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Cannot determine script directory.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
script_dir <- dirname(script_path)

required_packages <- c("dplyr", "tidyr", "rio")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R package(s): ", paste(missing_packages, collapse = ", "))
}

source(file.path(script_dir, "01_dhu_analysis.R"))
source(file.path(script_dir, "02_dhu_analysis_merge_rep.R"))

# These post-merge rules are the rules used to create the current HeLa
# 14_merge_DHU_site_rep tables. They intentionally remain explicit here so
# the published workflow reproduces the data tables without rendering QMDs.
assign_classic_sprinzl_pos <- function(position, sprinzl_position) {
  sp <- as.character(sprinzl_position)
  dplyr::case_when(
    sp == "20a" ~ "20a",
    sp %in% c("21", "22") ~ "20a",
    sp == "17" ~ "16",
    sp == "16" ~ "16",
    sp %in% c("19", "20") ~ "20",
    sp == "47" ~ "47",
    position == 17 ~ "16",
    position == 16 ~ "16",
    position == 19 ~ "20",
    position %in% c(21, 22) ~ "20a",
    position %in% c(46, 48) ~ "47",
    position == 47 ~ "47",
    position == 20 ~ "20",
    TRUE ~ NA_character_
  )
}

is_classic_sprinzl <- function(position, sprinzl_position) {
  !is.na(assign_classic_sprinzl_pos(position, sprinzl_position))
}

is_mt_trna_chr <- function(chr) grepl("mt-tRNA", chr, fixed = TRUE)
mt_single_rep_positions <- c(16, 17, 19, 20, 21, 46, 47, 48, 49)

is_mt_single_rep_exempt <- function(chr, position) {
  is_mt_trna_chr(chr) & position %in% mt_single_rep_positions
}

is_classic_for_filter <- function(chr, position, sprinzl_position) {
  is_classic_sprinzl(position, sprinzl_position) & !is_mt_trna_chr(chr)
}

get_rep_dhu_keys <- function(rep_file) {
  if (!file.exists(rep_file)) return(character())
  rep_df <- read.csv(rep_file)
  if (!all(c("chr", "position", "is_dhu_site") %in% names(rep_df))) {
    return(character())
  }
  rep_df |>
    dplyr::filter(is_dhu_site %in% TRUE) |>
    dplyr::mutate(key = paste(chr, position, sep = "|")) |>
    dplyr::pull(key) |>
    unique()
}

apply_mt_exempt_or_dhu_rule <- function(merged_data, rep1_file, rep2_file) {
  if (!all(c("chr", "position", "is_dhu_site") %in% names(merged_data))) {
    return(merged_data)
  }
  exempt_idx <- is_mt_single_rep_exempt(merged_data$chr, merged_data$position)
  if (!any(exempt_idx)) return(merged_data)
  row_keys <- paste(merged_data$chr, merged_data$position, sep = "|")
  rep_dhu_keys <- union(get_rep_dhu_keys(rep1_file), get_rep_dhu_keys(rep2_file))
  promote_idx <- exempt_idx & row_keys %in% rep_dhu_keys &
    !(merged_data$is_dhu_site %in% TRUE)
  if (any(promote_idx)) {
    merged_data$is_dhu_site[promote_idx] <- TRUE
    if ("site_name" %in% names(merged_data)) {
      merged_data$site_name[promote_idx] <- paste0(
        merged_data$chr[promote_idx], "-", merged_data$position[promote_idx]
      )
    }
  }
  merged_data
}

filter_oligo_dhu_chr <- function(df) {
  if (!"chr" %in% names(df)) return(df)
  df[!grepl("DHU", df$chr, fixed = TRUE), , drop = FALSE]
}

filter_nonclassic_both_rep <- function(df) {
  if (!all(c("sprinzl_position", "data_source") %in% names(df))) return(df)
  classic <- is_classic_for_filter(df$chr, df$position, df$sprinzl_position)
  drop_idx <- !classic & df$data_source != "both_replicates" &
    !is_mt_single_rep_exempt(df$chr, df$position)
  df[!drop_idx, , drop = FALSE]
}

filter_nonclassic_bg_mutation <- function(df, max_bg_mutation = 0.02) {
  req_cols <- c("mutation_rate", "correct_mutation", "sprinzl_position", "position")
  if (!all(req_cols %in% names(df)) || !"is_dhu_site" %in% names(df)) return(df)
  classic <- is_classic_for_filter(df$chr, df$position, df$sprinzl_position)
  bg_mutation <- df$mutation_rate - df$correct_mutation
  fail_idx <- !classic & df$is_dhu_site %in% TRUE &
    (is.na(bg_mutation) | bg_mutation > max_bg_mutation)
  if (any(fail_idx)) {
    df$is_dhu_site[fail_idx] <- FALSE
    if ("site_name" %in% names(df)) df$site_name[fail_idx] <- NA_character_
  }
  df
}

filter_nonclassic_not_dhu <- function(df) {
  if (!all(c("sprinzl_position", "position", "is_dhu_site") %in% names(df))) {
    return(df)
  }
  classic <- is_classic_for_filter(df$chr, df$position, df$sprinzl_position)
  df[!(!classic & !(df$is_dhu_site %in% TRUE)), , drop = FALSE]
}

filter_nonclassic_correct_mutation <- function(
  df,
  min_correct_mutation_cyto = 0.10,
  min_correct_mutation_mt = 0.10
) {
  req_cols <- c("correct_mutation", "sprinzl_position", "position", "chr")
  if (!all(req_cols %in% names(df))) return(df)
  classic <- is_classic_for_filter(df$chr, df$position, df$sprinzl_position)
  min_thr <- ifelse(
    is_mt_trna_chr(df$chr), min_correct_mutation_mt, min_correct_mutation_cyto
  )
  drop_idx <- !classic & (is.na(df$correct_mutation) | df$correct_mutation <= min_thr)
  df[!drop_idx, , drop = FALSE]
}

apply_post_merge_filters <- function(df) {
  df |>
    filter_oligo_dhu_chr() |>
    filter_nonclassic_both_rep() |>
    filter_nonclassic_bg_mutation() |>
    filter_nonclassic_not_dhu() |>
    filter_nonclassic_correct_mutation()
}

get_rep_venn_summary <- function(rep1_file, rep2_file) {
  get_sites <- function(path) {
    df <- read.csv(path)
    unique(paste(df$chr[df$is_dhu_site %in% TRUE],
                 df$position[df$is_dhu_site %in% TRUE], sep = "-"))
  }
  rep1_sites <- get_sites(rep1_file)
  rep2_sites <- get_sites(rep2_file)
  data.frame(
    rep1 = sub("\\.txt\\.csv$", "", basename(rep1_file)),
    rep2 = sub("\\.txt\\.csv$", "", basename(rep2_file)),
    rep1_only = length(setdiff(rep1_sites, rep2_sites)),
    rep2_only = length(setdiff(rep2_sites, rep1_sites)),
    both = length(intersect(rep1_sites, rep2_sites)),
    union = length(union(rep1_sites, rep2_sites)),
    stringsAsFactors = FALSE
  )
}

design <- read.csv(
  sample_sheet_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
required_columns <- c(
  "sample_id", "role", "background_set", "merge_group", "replicate",
  "summary_name", "summary_group"
)
missing_columns <- setdiff(required_columns, names(design))
if (length(missing_columns) > 0) {
  stop("Sample sheet is missing column(s): ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(design$sample_id)) stop("sample_id values must be unique.")
if (!all(design$role %in% c("background", "treatment"))) {
  stop("role must be either background or treatment.")
}

annotated_candidates <- list.files(
  annotated_dir,
  pattern = "\\.sprinzl\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)
resolve_annotated_file <- function(sample_id) {
  expected_name <- paste0(
    sample_id,
    "-merged-filtered-revcomp-filtered.R1-dedup-merged.sprinzl.txt"
  )
  matches <- annotated_candidates[basename(annotated_candidates) == expected_name]
  if (length(matches) != 1) {
    stop(
      "Expected exactly one annotated file named ", expected_name,
      " below ", annotated_dir, "; found ", length(matches), "."
    )
  }
  matches[[1]]
}
design$annotated_file <- vapply(
  design$sample_id,
  resolve_annotated_file,
  character(1)
)

background_rows <- design[design$role == "background", , drop = FALSE]
treatment_rows <- design[design$role == "treatment", , drop = FALSE]
if (nrow(background_rows) == 0) stop("No background rows in sample sheet.")
if (nrow(treatment_rows) == 0) stop("No treatment rows in sample sheet.")
if (any(is.na(treatment_rows$background_set))) {
  stop("Every treatment row must define background_set.")
}
if (any(is.na(treatment_rows$merge_group))) {
  stop("Every treatment row must define merge_group.")
}

dir.create(background_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(merged_output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Background correction: 11.anno_sprinzl -> 13_calu_bg_tRNA ===\n")
for (i in seq_len(nrow(treatment_rows))) {
  row <- treatment_rows[i, , drop = FALSE]
  background_files <- background_rows$annotated_file[
    background_rows$background_set == row$background_set
  ]
  if (length(background_files) == 0) {
    stop(
      "No background samples found for background_set: ",
      row$background_set
    )
  }

  cat("\nProcessing ", row$sample_id, " against ",
      length(background_files), " background sample(s)\n", sep = "")
  corrected <- process_dhu_with_background(
    main_file_path = row$annotated_file,
    background_file_path = background_files,
    min_mean_depth = 10,
    mutation_threshold = 0.02,
    stop_threshold = 0,
    filter_only_t = TRUE,
    change_low_cov_to_0 = FALSE,
    site_depth_to_0 = 10,
    save_data = FALSE,
    match_by_chr = TRUE,
    filter_mutation = TRUE,
    filter_stop = TRUE,
    filter_signal = TRUE,
    filter_t2c = TRUE,
    signal_threshold = 0.03,
    t2c_threshold = 0.3,
    consecutive_t_correction = TRUE,
    consecutive_t_mut_threshold = 0.02,
    clear_all_if_filter = FALSE,
    filter_dhu_site = TRUE,
    chr_filter_pattern = NULL,
    max_chr_length = 200,
    dhu_site_min_depth = 50
  )
  # Preserve the original QMD workflow's sample_name convention exactly.
  corrected$sample_name <- paste0(row$sample_id, ".txt")
  write.csv(
    corrected,
    file.path(background_output_dir, paste0(row$sample_id, ".txt.csv")),
    row.names = FALSE
  )
}

cat("\n=== Replicate merge: 13_calu_bg_tRNA -> 14_merge_DHU_site_rep ===\n")
merge_groups <- unique(treatment_rows$merge_group)
venn_summaries <- list()
for (group_name in merge_groups) {
  group_rows <- treatment_rows[
    treatment_rows$merge_group == group_name,
    ,
    drop = FALSE
  ]
  group_rows <- group_rows[order(group_rows$replicate), , drop = FALSE]
  if (!nrow(group_rows) %in% c(1, 2)) {
    stop(
      "merge_group '", group_name,
      "' must contain one sample or exactly two replicates."
    )
  }

  rep1_file <- file.path(
    background_output_dir,
    paste0(group_rows$sample_id[[1]], ".txt.csv")
  )
  output_file <- file.path(merged_output_dir, paste0(group_name, ".txt.csv"))

  if (nrow(group_rows) == 1) {
    single_data <- rio::import(rep1_file)
    single_data <- apply_post_merge_filters(single_data)
    rio::export(single_data, output_file)
    next
  }
  if (!identical(as.integer(group_rows$replicate), c(1L, 2L))) {
    stop("Two-sample merge_group '", group_name, "' must contain replicate 1 and 2.")
  }
  rep2_file <- file.path(
    background_output_dir,
    paste0(group_rows$sample_id[[2]], ".txt.csv")
  )

  merged_data <- merge_dhu_replicates(
    replicate1_file = rep1_file,
    replicate2_file = rep2_file,
    keep_in_one_file = TRUE,
    output_file = NULL,
    dhu_site_rule = "and",
    merged_sample_name = paste0(group_name, ".txt"),
    verbose = TRUE
  )
  merged_data <- apply_mt_exempt_or_dhu_rule(merged_data, rep1_file, rep2_file)
  merged_data <- apply_post_merge_filters(merged_data)
  rio::export(merged_data, output_file)

  # This helper table is present in the published HeLa folder 14.
  if (identical(group_name, "HeLa_tRNA_IP_Induro_Mg")) {
    met_data <- merged_data[
      grepl("Met-CAT", merged_data$chr, fixed = TRUE),
      ,
      drop = FALSE
    ]
    rio::export(
      met_data,
      file.path(merged_output_dir, paste0(group_name, "_met_chr.csv"))
    )
  }

  summary_name <- unique(group_rows$summary_name)
  if (length(summary_name) != 1 || is.na(summary_name)) {
    stop("Each merge_group must have one non-empty summary_name.")
  }
  summary_row <- get_rep_venn_summary(rep1_file, rep2_file)
  summary_group <- unique(group_rows$summary_group)
  if (length(summary_group) != 1 || is.na(summary_group)) {
    stop("Each two-replicate merge_group must have one non-empty summary_group.")
  }
  summary_row$rt_condition <- summary_group
  summary_row <- summary_row[, c(
    "rt_condition", "rep1", "rep2", "rep1_only", "rep2_only", "both", "union"
  )]
  venn_summaries[[summary_name]] <- rbind(
    venn_summaries[[summary_name]], summary_row
  )
}

for (summary_name in names(venn_summaries)) {
  write.csv(
    venn_summaries[[summary_name]],
    file.path(merged_output_dir, paste0(summary_name, "_venn_summary.csv")),
    row.names = FALSE
  )
}

cat("\nPost-processing completed. Final files: ", merged_output_dir, "\n", sep = "")
