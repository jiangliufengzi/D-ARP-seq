# DHU replicate merging - step 2
# Author: generated for the DHU analysis workflow
#
# Functions for merging two replicate files produced by dhu_analysis_step1.R
# Main function: merge_dhu_replicates() - merge two background-corrected DHU replicates

#' Merge two background-corrected DHU replicate files
#'
#' Merge two replicate tables that have already been processed and
#' background-corrected by 01_dhu_analysis.R.
#' For positions present in both replicates, numeric signal columns are
#' averaged; structural columns (chr_length and similar) keep the rep1 value;
#' is_dhu_site is combined with dhu_site_rule; site_name is rebuilt from the
#' merged is_dhu_site; sample_name is set to a single merged name.
#' Positions present in only one replicate are controlled by keep_in_one_file.
#'
#' @param replicate1_file Path to the first background-corrected replicate file
#' @param replicate2_file Path to the second background-corrected replicate file
#' @param keep_in_one_file Logical. If TRUE, keep positions found in only one
#'   replicate using their original values; if FALSE, drop them (default: TRUE)
#' @param output_file Optional path for the merged table (default: NULL, do not save)
#' @param replicate1_name Label for replicate 1 (default: "rep1")
#' @param replicate2_name Label for replicate 2 (default: "rep2")
#' @param dhu_site_rule Merge rule for is_dhu_site: "and" requires TRUE in both
#'   replicates (default, strictest / most reproducible), "or" accepts TRUE in
#'   either replicate, "rep1" keeps the rep1 call. Single-replicate positions
#'   are always FALSE under "and" (reproducibility cannot be assessed)
#' @param merged_sample_name Value for the merged sample_name (default NULL:
#'   inferred from rep1 sample_name after stripping rep1/rep2 tags)
#' @param take_first_cols Numeric columns that should be identical across
#'   replicates and therefore take the rep1 value instead of averaging
#'   (default c("chr_length","seq_length"))
#' @param verbose Logical. If TRUE, print progress (default: TRUE)
#' @return Merged replicate data frame
#' @export
merge_dhu_replicates <- function(replicate1_file,
                                replicate2_file,
                                keep_in_one_file = TRUE,
                                output_file = NULL,
                                replicate1_name = "rep1",
                                replicate2_name = "rep2",
                                dhu_site_rule = c("and", "or", "rep1"),
                                merged_sample_name = NULL,
                                take_first_cols = c("chr_length", "seq_length", "sprinzl_position", "ref_base", "strand"),
                                verbose = TRUE) {
  dhu_site_rule <- match.arg(dhu_site_rule)
  
  # Load required R packages
  require(dplyr)
  require(tidyr)
  require(rio)
  
  if(verbose) {
    cat("=== DHU Replicate Merger (Step 2) ===\n")
    cat("Replicate 1:", replicate1_file, "\n")
    cat("Replicate 2:", replicate2_file, "\n")
    cat("Keep single-replicate positions:", keep_in_one_file, "\n")
  }
  
  # Load replicate files
  if(verbose) cat("Loading replicate files...\n")
  rep1_data <- normalize_replicate_import(rio::import(replicate1_file))
  rep2_data <- normalize_replicate_import(rio::import(replicate2_file))
  
  if(verbose) {
    cat("Replicate 1 data:", nrow(rep1_data), "rows,", ncol(rep1_data), "columns\n")
    cat("Replicate 2 data:", nrow(rep2_data), "rows,", ncol(rep2_data), "columns\n")
  }
  
  # Validate the input tables
  validate_replicate_data(rep1_data, rep2_data, verbose)
  
  # Tag each table with its replicate source
  rep1_data$replicate_source <- replicate1_name
  rep2_data$replicate_source <- replicate2_name
  
  # Join keys (chr + position)
  join_keys <- c("chr", "position")

  # These columns are handled separately (not averaged / not generic rep1-first)
  special_cols <- intersect(c("is_dhu_site", "site_name", "sample_name"), names(rep1_data))

  # Numeric columns to average (exclude keys, source, structural, and special columns)
  numeric_cols <- identify_numeric_columns(
    rep1_data,
    exclude_cols = c(join_keys, "replicate_source", take_first_cols, special_cols)
  )

  if(verbose) {
    cat("Joining by:", paste(join_keys, collapse = ", "), "\n")
    cat("Numeric columns to average:", length(numeric_cols), "columns\n")
    cat("is_dhu_site merge rule:", dhu_site_rule, "\n")
  }

  # Full join so positions unique to either replicate are retained
  merged_data <- full_join(rep1_data, rep2_data,
                          by = join_keys,
                          suffix = c(".rep1", ".rep2"))

  # Record which replicate(s) contributed each position
  merged_data <- merged_data %>%
    mutate(
      # Assign data_source
      data_source = case_when(
        !is.na(replicate_source.rep1) & !is.na(replicate_source.rep2) ~ "both_replicates",
        !is.na(replicate_source.rep1) & is.na(replicate_source.rep2) ~ replicate1_name,
        is.na(replicate_source.rep1) & !is.na(replicate_source.rep2) ~ replicate2_name,
        TRUE ~ "unknown"
      )
    )

  # (1) Numeric signal columns: average when both replicates have a value
  if(verbose) cat("Calculating averaged values...\n")

  for(col in numeric_cols) {
    col_rep1 <- paste0(col, ".rep1")
    col_rep2 <- paste0(col, ".rep2")

    merged_data[[col]] <- calculate_merged_value(
      merged_data[[col_rep1]],
      merged_data[[col_rep2]],
      keep_single = keep_in_one_file
    )

    merged_data[[col_rep1]] <- NULL
    merged_data[[col_rep2]] <- NULL
  }

  # (2) Structural columns (should match): take rep1, then rep2; do not average
  for(col in intersect(take_first_cols, names(rep1_data))) {
    col_rep1 <- paste0(col, ".rep1")
    col_rep2 <- paste0(col, ".rep2")
    if(col_rep1 %in% names(merged_data) && col_rep2 %in% names(merged_data)) {
      merged_data[[col]] <- coalesce(merged_data[[col_rep1]], merged_data[[col_rep2]])
      merged_data[[col_rep1]] <- NULL
      merged_data[[col_rep2]] <- NULL
    }
  }

  # (3) Remaining metadata columns (ref_base/strand/sprinzl_position/*_yes_or_no): rep1 first
  non_numeric_cols <- setdiff(
    names(rep1_data),
    c(join_keys, numeric_cols, take_first_cols, special_cols, "replicate_source")
  )
  for(col in non_numeric_cols) {
    col_rep1 <- paste0(col, ".rep1")
    col_rep2 <- paste0(col, ".rep2")

    if(col_rep1 %in% names(merged_data) && col_rep2 %in% names(merged_data)) {
      merged_data[[col]] <- ifelse(!is.na(merged_data[[col_rep1]]),
                                  merged_data[[col_rep1]],
                                  merged_data[[col_rep2]])
      merged_data[[col_rep1]] <- NULL
      merged_data[[col_rep2]] <- NULL
    }
  }

  # (4) Merge is_dhu_site by rule (missing values in one replicate count as FALSE)
  if(all(c("is_dhu_site.rep1", "is_dhu_site.rep2") %in% names(merged_data))) {
    f1 <- coalesce(merged_data[["is_dhu_site.rep1"]], FALSE)
    f2 <- coalesce(merged_data[["is_dhu_site.rep2"]], FALSE)
    merged_data$is_dhu_site <- switch(
      dhu_site_rule,
      and  = f1 & f2,
      or   = f1 | f2,
      rep1 = coalesce(merged_data[["is_dhu_site.rep1"]], merged_data[["is_dhu_site.rep2"]])
    )
    merged_data[["is_dhu_site.rep1"]] <- NULL
    merged_data[["is_dhu_site.rep2"]] <- NULL
  }

  # (5) Rebuild site_name from the merged is_dhu_site so the two stay consistent
  for(sn in c("site_name.rep1", "site_name.rep2")) {
    if(sn %in% names(merged_data)) merged_data[[sn]] <- NULL
  }
  if("is_dhu_site" %in% names(merged_data)) {
    merged_data$site_name <- ifelse(
      merged_data$is_dhu_site,
      paste0(merged_data$chr, "-", merged_data$position),
      NA_character_
    )
  }

  # (6) Collapse sample_name to a single merged value
  for(sn in c("sample_name.rep1", "sample_name.rep2")) {
    if(sn %in% names(merged_data)) merged_data[[sn]] <- NULL
  }
  if("sample_name" %in% names(rep1_data)) {
    if(is.null(merged_sample_name)) {
      s1 <- as.character(rep1_data$sample_name[1])
      # Strip a delimiter-prefixed rep1/rep2 tag and anything after it.
      # Require '.', '_', or '-' before "rep" so names that contain "rep" mid-string
      # are not truncated.
      merged_sample_name <- sub("[._-]rep[12].*$", "", s1)
      if(is.na(merged_sample_name) || merged_sample_name == "") merged_sample_name <- s1
    }
    merged_data$sample_name <- merged_sample_name
  }

  # Drop replicate-source columns
  merged_data$replicate_source.rep1 <- NULL
  merged_data$replicate_source.rep2 <- NULL
  
  # Optionally drop positions present in only one replicate
  if(!keep_in_one_file) {
    original_rows <- nrow(merged_data)
    merged_data <- merged_data %>%
      filter(data_source == "both_replicates")
    if(verbose) {
      cat("Removed", original_rows - nrow(merged_data), "positions present in only one replicate\n")
    }
  }

  # Column order: chr, position, sprinzl_position first; remaining columns follow rep1
  if("sprinzl_position" %in% names(merged_data)) {
    lead_cols <- intersect(c("chr", "position", "sprinzl_position"), names(merged_data))
    tail_cols <- setdiff(names(merged_data), lead_cols)
    ordered_tail <- c(
      intersect(names(rep1_data), tail_cols),
      setdiff(tail_cols, names(rep1_data))
    )
    merged_data <- merged_data[, c(lead_cols, ordered_tail), drop = FALSE]
  }
  
  # Summarize merge counts
  summary_stats <- merged_data %>%
    group_by(data_source) %>%
    summarise(count = n(), .groups = 'drop')
  
  if(verbose) {
    cat("\nMerge summary:\n")
    print(summary_stats)
    cat("Final merged data:", nrow(merged_data), "rows,", ncol(merged_data), "columns\n")
  }
  
  # Save when an output path is provided
  if(!is.null(output_file)) {
    if(verbose) cat("Saving merged data to:", output_file, "\n")
    rio::export(merged_data, output_file)
  }
  
  if(verbose) cat("DHU replicate merging completed successfully!\n")
  
  return(merged_data)
}

#' Validate replicate table structure and compatibility
#'
#' @param rep1_data First replicate data frame
#' @param rep2_data Second replicate data frame
#' @param verbose Whether to print progress
validate_replicate_data <- function(rep1_data, rep2_data, verbose = TRUE) {
  
  # Check required columns
  required_cols <- c("chr", "position")
  
  missing_rep1 <- setdiff(required_cols, names(rep1_data))
  missing_rep2 <- setdiff(required_cols, names(rep2_data))
  
  if(length(missing_rep1) > 0) {
    stop("Replicate 1 missing required columns: ", paste(missing_rep1, collapse = ", "))
  }
  
  if(length(missing_rep2) > 0) {
    stop("Replicate 2 missing required columns: ", paste(missing_rep2, collapse = ", "))
  }
  
  # Check column overlap
  common_cols <- intersect(names(rep1_data), names(rep2_data))
  if(verbose) {
    cat("Common columns between replicates:", length(common_cols), "\n")
  }
  
  # Report columns unique to one replicate
  rep1_only <- setdiff(names(rep1_data), names(rep2_data))
  rep2_only <- setdiff(names(rep2_data), names(rep1_data))
  
  if(length(rep1_only) > 0 && verbose) {
    cat("Columns only in replicate 1:", paste(rep1_only, collapse = ", "), "\n")
  }
  
  if(length(rep2_only) > 0 && verbose) {
    cat("Columns only in replicate 2:", paste(rep2_only, collapse = ", "), "\n")
  }
}

#' Identify numeric columns that should be averaged
#'
#' @param data Data frame to inspect
#' @param exclude_cols Columns to leave out of averaging
#' @return Character vector of numeric column names
normalize_replicate_import <- function(data) {
  char_cols <- c(
    "chr", "sprinzl_position", "ref_base", "strand",
    "mutation_yes_or_no", "stop_yes_or_no", "site_name", "sample_name"
  )
  for (col in intersect(char_cols, names(data))) {
    data[[col]] <- as.character(data[[col]])
  }
  if ("is_dhu_site" %in% names(data)) {
    data$is_dhu_site <- as.logical(data$is_dhu_site)
  }
  data
}

identify_numeric_columns <- function(data, exclude_cols = c()) {
  
  # All numeric columns
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  
  # Drop excluded columns
  numeric_cols <- setdiff(numeric_cols, exclude_cols)
  
  return(numeric_cols)
}

#' Combine a pair of replicate-column values
#'
#' @param val1 Values from replicate 1
#' @param val2 Values from replicate 2
#' @param keep_single Whether to keep values present in only one replicate
#' @return Combined value vector
calculate_merged_value <- function(val1, val2, keep_single = TRUE) {
  
  # Initialize the result vector
  result <- rep(NA, length(val1))
  
  # Both values present: take the mean
  both_present <- !is.na(val1) & !is.na(val2)
  result[both_present] <- (val1[both_present] + val2[both_present]) / 2
  
  if(keep_single) {
    # Only val1 is present
    only_val1 <- !is.na(val1) & is.na(val2)
    result[only_val1] <- val1[only_val1]
    
    # Only val2 is present
    only_val2 <- is.na(val1) & !is.na(val2)
    result[only_val2] <- val2[only_val2]
  }
  
  return(result)
}

#' Merge multiple replicate pairs in batch
#'
#' @param replicate_pairs List of lists, each with replicate1_file and replicate2_file
#' @param output_dir Optional directory for merged files
#' @param keep_in_one_file Passed through to merge_dhu_replicates
#' @param verbose Whether to print progress
#' @return List of merged data frames
batch_merge_replicates <- function(replicate_pairs,
                                  output_dir = NULL,
                                  keep_in_one_file = TRUE,
                                  verbose = TRUE,
                                  ...) {
  
  if(verbose) {
    cat("=== Batch DHU Replicate Merger ===\n")
    cat("Processing", length(replicate_pairs), "replicate pairs\n")
  }
  
  merged_results <- list()
  
  for(i in seq_along(replicate_pairs)) {
    pair <- replicate_pairs[[i]]
    
    if(verbose) cat("\nProcessing pair", i, "of", length(replicate_pairs), "\n")
    
    # Build an output path when output_dir is set
    output_file <- NULL
    if(!is.null(output_dir)) {
      base_name <- paste0("merged_pair_", i, "_", Sys.Date(), ".csv")
      output_file <- file.path(output_dir, base_name)
    }
    
    # Merge the current replicate pair
    merged_data <- merge_dhu_replicates(
      replicate1_file = pair$replicate1_file,
      replicate2_file = pair$replicate2_file, 
      keep_in_one_file = keep_in_one_file,
      output_file = output_file,
      replicate1_name = paste0("rep1_pair", i),
      replicate2_name = paste0("rep2_pair", i),
      verbose = verbose,
      ...
    )
    
    merged_results[[i]] <- merged_data
  }
  
  if(verbose) cat("\nBatch merging completed for", length(replicate_pairs), "pairs!\n")
  
  return(merged_results)
}
