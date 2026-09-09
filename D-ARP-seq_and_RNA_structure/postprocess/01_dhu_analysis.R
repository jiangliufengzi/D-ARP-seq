# DHU Analysis Functions
# Author: Generated for DHU figure1 analysis
# 
# Contains three main functions:
# 1. process_dhu_data() - Single file DHU analysis
# 2. process_dhu_with_background() - DHU analysis with background correction
# 3. batch_dhu_with_background() - Batch DHU analysis with shared background

dhu_safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  max(x)
}

dhu_prepare_background_for_join <- function(bg_data_list, join_keys) {
  if (length(bg_data_list) == 0) {
    stop("No background data available for joining.")
  }
  
  bg_data <- dplyr::bind_rows(bg_data_list)
  bg_data <- dplyr::group_by(bg_data, dplyr::across(dplyr::all_of(join_keys)))
  bg_data <- dplyr::summarise(
    bg_data,
    bg_mutation_rate = mean(bg_mutation_rate, na.rm = TRUE),
    bg_normalized_mutation_rate = mean(bg_normalized_mutation_rate, na.rm = TRUE),
    bg_stop_rate = mean(bg_stop_rate, na.rm = TRUE),
    bg_normalized_stop_rate = mean(bg_normalized_stop_rate, na.rm = TRUE),
    bg_signal_value = mean(bg_signal_value, na.rm = TRUE),
    .groups = "drop"
  )
  dplyr::mutate(
    bg_data,
    dplyr::across(
      dplyr::starts_with("bg_"),
      ~ ifelse(is.nan(.x), NA_real_, .x)
    )
  )
}

dhu_apply_dhu_site_flag <- function(data,
                                    filter_mutation = FALSE,
                                    filter_stop = FALSE,
                                    filter_signal = FALSE,
                                    filter_t2c = FALSE,
                                    mutation_threshold = 0.05,
                                    stop_threshold = 0.05,
                                    signal_threshold = 0.05,
                                    t2c_threshold = 0.66,
                                    dhu_site_min_depth = 50,
                                    chr_filter_pattern = NULL,
                                    hard_filter = FALSE) {
  data <- dplyr::mutate(
    data,
    is_dhu_site = ref_base == "T" &
      depth >= dhu_site_min_depth &
      position > 10 &
      correct_mutation > 0 &
      correct_signal > 0 &
      correct_signal < 1 &
      (!filter_t2c | t2c_dominance >= t2c_threshold) &
      (!filter_mutation | normalized_correct_mutation >= mutation_threshold) &
      (!filter_stop | normalized_correct_stop >= stop_threshold) &
      (!filter_signal | correct_signal >= signal_threshold)
  )
  
  if (!is.null(chr_filter_pattern) && chr_filter_pattern != "") {
    data <- dplyr::mutate(data, is_dhu_site = is_dhu_site & grepl(chr_filter_pattern, chr))
  }
  
  data <- dplyr::mutate(data, site_name = ifelse(is_dhu_site, paste0(chr, "-", position), NA_character_))
  
  if (hard_filter) {
    data <- dplyr::filter(data, is_dhu_site)
  }
  
  data
}

dhu_zero_first_positions <- function(data, clear_first_positions = 0) {
  if (is.null(clear_first_positions) || clear_first_positions <= 0 || !"position" %in% names(data)) {
    return(data)
  }
  
  rows <- which(data$position <= clear_first_positions)
  if (length(rows) == 0) {
    return(data)
  }
  
  cols_to_zero <- intersect(
    c(
      "correct_mutation", "correct_normalized_mutation",
      "correct_stop", "correct_normalized_stop",
      "correct_signal",
      "normalized_correct_mutation", "normalized_correct_stop",
      "normalized_correct_signal"
    ),
    names(data)
  )
  
  for (col in cols_to_zero) {
    data[[col]][rows] <- 0
  }
  
  if ("is_dhu_site" %in% names(data)) {
    data$is_dhu_site[rows] <- FALSE
  }
  if ("site_name" %in% names(data)) {
    data$site_name[rows] <- NA_character_
  }
  
  data
}

#' Process DHU sequencing data to calculate normalized mutation and stop rates
#'
#' @param file_path Path to the merged data file
#' @param min_mean_depth Minimum mean depth threshold for chromosome filtering (default: 10)
#' @param mutation_threshold Threshold for mutation signal detection (default: 0.05)
#' @param stop_threshold Threshold for stop signal detection (default: 0.05)
#' @param filter_only_t Whether to filter only T reference bases (default: FALSE). If TRUE, sets mutation and stop rates to 0 for non-T reference bases
#' @param change_low_cov_to_0 Whether to apply low coverage filtering (default: FALSE). If TRUE, sets depth<site_depth_to_0 to 0 and position 1-3 stop rates to 0
#' @param site_depth_to_0 Depth threshold for low coverage filtering (default: 50). Sites with depth below this value are set to 0 when change_low_cov_to_0 is TRUE
#' @param save_data Whether to save processed data (default: FALSE)
#' @param data_output_file Output file path for processed data (optional)
#' @param filter_mutation Whether to filter positions with low normalized mutation rates (default: FALSE). If TRUE, positions with normalized_mutation_rate < mutation_threshold will have signal data set to 0
#' @param filter_stop Whether to filter positions with low normalized stop rates (default: FALSE). If TRUE, positions with normalized_stop_rate < stop_threshold will have signal data set to 0
#' @param filter_signal Whether to filter positions with low signal values (default: FALSE). If TRUE, positions with signal_value < signal_threshold will have signal data set to 0
#' @param signal_threshold Threshold for signal filtering (default: 0.05). Used when filter_signal is TRUE
#' @param clear_all_if_filter Whether to clear all data including technical columns when filtering (default: FALSE). If TRUE, clears all data; if FALSE, only clears main signal columns
#' @param filter_dhu_site Whether to apply final DHU site filtering (default: FALSE). If TRUE, applies filter_dhu conditions: ref_base=='T', depth>=dhu_site_min_depth, normalized_mutation_rate>=0.05, signal_value>=0.08, position>10, position<70, signal_value<1, tRNA chromosomes only
#' @param chr_filter_pattern Pattern to filter chromosome names (default: "Homo_sapiens_tRNA"). Used when filter_dhu_site=TRUE. Set to NULL or empty string ("") to skip chromosome filtering
#' @param dhu_site_min_depth Minimum depth required for final DHU site filtering (default: 50)
#' @return Processed data frame with normalized rates and classifications
process_dhu_data <- function(file_path, 
                             min_mean_depth = 10, 
                             mutation_threshold = 0.05, 
                             stop_threshold = 0.05,
                             filter_only_t = FALSE,
                             change_low_cov_to_0 = FALSE,
                             site_depth_to_0 = 50,
                             save_data = FALSE,
                             data_output_file = NULL,
                             filter_mutation = FALSE,
                             filter_stop = FALSE,
                             filter_signal = FALSE,
                             signal_threshold = 0.05,
                             clear_all_if_filter = FALSE,
                             filter_dhu_site = FALSE,
                             chr_filter_pattern = NULL,
                             dhu_site_min_depth = 50) {
  
  # Load required libraries
  require(dplyr)
  require(tidyr)
  require(rio)
  
  # Import data
  df <- rio::import(file_path)
  
  ########################################################################
  # Filter out low coverage chromosomes
  ########################################################################
  low_cov <- df %>%
    group_by(chr) %>%
    summarise(mean_depth = mean(depth), .groups = 'drop') %>%
    filter(mean_depth < min_mean_depth)
  
  df <- df %>%
    filter(!chr %in% low_cov$chr)
  
  ########################################################################
  # Fill missing positions for each chromosome
  ########################################################################
  fill_missing_positions <- function(df) {
    has_sprinzl_position <- "sprinzl_position" %in% names(df)
    
    # Get unique chromosome information
    chr_info <- df %>%
      select(chr, chr_length, strand) %>%
      distinct()
    
    # Generate complete position sequences for each chr-strand combination
    complete_positions <- chr_info %>%
      rowwise() %>%
      do({
        data.frame(
          chr = .$chr,
          position = 1:.$chr_length,
          strand = .$strand,
          chr_length = .$chr_length
        )
      }) %>%
      ungroup()
    
    # Get default values for missing positions
    default_values <- data.frame(
      ref_base = "N",  # Default for missing positions
      stop_value = 0.0,
      depth = 0,
      A = 0, C = 0, G = 0, `T` = 0, N = 0,
      P = 0, PM = 0, S = 0, SM = 0,
      mutation_rate = 0.0,
      signal_value = 0.0  # Add default value for signal_value
    )
    
    # Combine complete positions with default values
    complete_data <- complete_positions %>%
      crossing(default_values)
    
    # Merge with original data, keeping original data priority
    output_cols <- c(
      "chr", "position",
      if (has_sprinzl_position) "sprinzl_position",
      "ref_base", "strand", "stop_value", "depth", "A", "C", "G", "T",
      "N", "P", "PM", "S", "SM", "mutation_rate", "signal_value", "chr_length"
    )
    
    result <- complete_data %>%
      left_join(df, by = c("chr", "position", "strand"), suffix = c("", ".orig")) %>%
      mutate(
        ref_base = coalesce(ref_base.orig, ref_base),
        stop_value = coalesce(stop_value.orig, stop_value),
        depth = coalesce(depth.orig, depth),
        A = coalesce(A.orig, A),
        C = coalesce(C.orig, C),
        G = coalesce(G.orig, G),
        `T` = coalesce(`T.orig`, `T`),
        N = coalesce(N.orig, N),
        P = coalesce(P.orig, P),
        PM = coalesce(PM.orig, PM),
        S = coalesce(S.orig, S),
        SM = coalesce(SM.orig, SM),
        mutation_rate = coalesce(mutation_rate.orig, mutation_rate),
        signal_value = coalesce(signal_value.orig, signal_value),
        chr_length = coalesce(chr_length.orig, chr_length)
      ) %>%
      select(dplyr::all_of(output_cols)) %>%
      arrange(chr, strand, position)
    
    return(result)
  }
  
  # Apply position filling
  df <- fill_missing_positions(df)
  
  ########################################################################
  # change_low_cov_site_depth_to_0
  ########################################################################
  if(change_low_cov_to_0) {
    df <- df %>%
      mutate(
        # Store original depth for consistent comparison
        original_depth = depth,
        # Apply filtering based on original depth values
        stop_value = ifelse(original_depth < site_depth_to_0, 0, stop_value),
        depth = ifelse(original_depth < site_depth_to_0, 0, depth),
        A = ifelse(original_depth < site_depth_to_0, 0, A),
        C = ifelse(original_depth < site_depth_to_0, 0, C),
        G = ifelse(original_depth < site_depth_to_0, 0, G),
        `T` = ifelse(original_depth < site_depth_to_0, 0, `T`),  # Use backticks to avoid confusion with TRUE
        N = ifelse(original_depth < site_depth_to_0, 0, N),
        P = ifelse(original_depth < site_depth_to_0, 0, P),
        PM = ifelse(original_depth < site_depth_to_0, 0, PM),
        S = ifelse(original_depth < site_depth_to_0, 0, S),
        SM = ifelse(original_depth < site_depth_to_0, 0, SM),
        mutation_rate = ifelse(original_depth < site_depth_to_0, 0, mutation_rate),
        signal_value = ifelse(original_depth < site_depth_to_0, 0, signal_value)  # Also set signal_value to 0 for low coverage
      ) %>%
      # Remove temporary column
      select(-original_depth)
  }
  
  ########################################################################
  # Recalculate mutation rate
  ########################################################################
  df <- df %>%
    mutate(
      # Determine reference base count based on ref_base
      ref_count = case_when(
        ref_base == "A" ~ A,
        ref_base == "C" ~ C,
        ref_base == "G" ~ G,
        ref_base == "T" ~ `T`,
        TRUE ~ 0  # Default for other cases like "N"
      ),
      # Calculate new mutation rate
      mutation_rate = case_when(
        depth == 0 ~ 0,  # If depth is 0, mutation rate is 0
        (depth - N) <= 0 ~ 0,  # If denominator is 0 or negative, mutation rate is 0
        
        TRUE ~ (depth - N - ref_count) / (depth - N)
      ),
      # Ensure mutation rate is between 0 and 1
      mutation_rate = pmax(0, pmin(1, mutation_rate))
    )
  
  ########################################################################
  # Calculate normalized mutation rate
  ########################################################################
  df <- df %>%
    group_by(chr) %>%
    mutate(mean_mutation_rate = mean(mutation_rate, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(normalized_mutation_rate = pmax(mutation_rate - mean_mutation_rate, 0)) %>%
    mutate(mutation_yes_or_no = ifelse(normalized_mutation_rate > mutation_threshold, "yes", "no"))
  
  ########################################################################
  # Calculate stop rate
  ########################################################################
  df <- df %>%
    group_by(chr, strand) %>%
    arrange(chr, strand, position) %>%
    mutate(
      # Get previous and next depth values
      prev_depth = lag(depth, default = 0),
      next_depth = lead(depth, default = 0),
      # Get sequence length for last 20 positions
      seq_length = max(position),
      # Calculate stop rate
      stop_rate = case_when(
        # position=1,2,3 set to 0 (controlled by change_low_cov_to_0)
        change_low_cov_to_0 & position == 1 ~ 0,
        change_low_cov_to_0 & position == 2 ~ 0,
        change_low_cov_to_0 & position == 3 ~ 0,
        # Last 20 positions set to 0
        position > (seq_length - 20) ~ 0,
        # When ref_base is T
        ref_base == "T" ~ ifelse(depth > 0 & next_depth > 0, (next_depth - prev_depth) / next_depth, 0),
        # when prev_depth is 0, keep stop stop_rate 0
        prev_depth == 0 ~ 0,
        # When ref_base is not T
        TRUE ~ ifelse(depth > 0, (depth - prev_depth) / depth, 0)
      )
    ) %>%
    # Ensure stop rate is non-negative
    mutate(stop_rate = pmax(stop_rate, 0)) %>%
    ungroup()
  
  ########################################################################
  # Calculate normalized stop rate
  ########################################################################
  df <- df %>%
    group_by(chr) %>%
    mutate(mean_stop_rate = mean(stop_rate, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(normalized_stop_rate = pmax(stop_rate - mean_stop_rate, 0)) %>%
    mutate(stop_yes_or_no = ifelse(normalized_stop_rate > stop_threshold, "yes", "no"))
  
  # Set first 3 positions of each chromosome to 0 for stop and mutation values
  # to exclude technical bias at chromosome starts (controlled by change_low_cov_to_0)
  if(change_low_cov_to_0) {
    df <- df %>%
      mutate(
        normalized_mutation_rate = ifelse(position <= 3, 0, normalized_mutation_rate),
        mutation_yes_or_no = ifelse(position <= 3, "no", mutation_yes_or_no)
      )
  }
  
  ########################################################################
  # Recalculate signal value
  ########################################################################
  df <- df %>%
    group_by(chr, strand) %>%
    arrange(chr, strand, position) %>%
    mutate(
      # Get previous and next depth values
      prev_S = lag(S, default = 0),
      prev_SM = lag(SM, default = 0),
      next_S = lead(S, default = 0),
      next_SM = lead(SM, default = 0),
      # Get sequence length for last 20 positions
      seq_length = max(position),
      # Calculate stop rate
      signal_value = case_when(
        # position=1,2,3 set to 0 (controlled by change_low_cov_to_0)
        change_low_cov_to_0 & position == 1 ~ 0,
        change_low_cov_to_0 & position == 2 ~ 0,
        change_low_cov_to_0 & position == 3 ~ 0,
        # Last 20 positions set to 0
        position > (seq_length - 20) ~ 0,
        # When ref_base is T
        ref_base == "T" ~ ifelse(
          depth > 0 & (depth + next_S + next_SM) > 0, 
          ((depth - P) + next_S + next_SM) / (depth + next_S + next_SM), 
          0
        ),
        # When ref_base is not T
        TRUE ~ ifelse(depth > 0, (depth - P) / depth, 0)
      )
    ) %>%
    # Ensure stop rate is non-negative
    mutate(signal_value = pmax(signal_value, 0)) %>%
    ungroup()
  
  ########################################################################
  # filter_only_t
  ########################################################################
  if(filter_only_t) {
    df <- df %>%
      mutate(
        # Set mutation-related values to 0 for non-T reference bases
        mutation_rate = ifelse(ref_base != "T", 0, mutation_rate),
        normalized_mutation_rate = ifelse(ref_base != "T", 0, normalized_mutation_rate),
        mutation_yes_or_no = ifelse(ref_base != "T", "no", mutation_yes_or_no),
        # Set stop-related values to 0 for non-T reference bases
        stop_rate = ifelse(ref_base != "T", 0, stop_rate),
        normalized_stop_rate = ifelse(ref_base != "T", 0, normalized_stop_rate),
        stop_yes_or_no = ifelse(ref_base != "T", "no", stop_yes_or_no),
        # Set signal_value and other columns to 0 for non-T reference bases
        signal_value = ifelse(ref_base != "T", 0, signal_value),
        P = ifelse(ref_base != "T", 0, P),
        PM = ifelse(ref_base != "T", 0, PM),
        S = ifelse(ref_base != "T", 0, S),
        SM = ifelse(ref_base != "T", 0, SM)
      )
  }

  ########################################################################
  # Apply mutation filtering if requested
  ########################################################################
  if(filter_mutation) {
    df <- df %>%
      mutate(
        # Store original mutation condition to avoid circular dependency
        mutation_filter_condition = normalized_mutation_rate < mutation_threshold,
        # Always clear main signal columns if normalized_mutation_rate < mutation_threshold
        mutation_rate = ifelse(mutation_filter_condition, 0, mutation_rate),
        normalized_mutation_rate = ifelse(mutation_filter_condition, 0, normalized_mutation_rate),
        stop_rate = ifelse(mutation_filter_condition, 0, stop_rate),
        normalized_stop_rate = ifelse(mutation_filter_condition, 0, normalized_stop_rate),
        signal_value = ifelse(mutation_filter_condition, 0, signal_value),
        # Only clear technical columns if clear_all_if_filter = TRUE
        stop_value = ifelse(clear_all_if_filter & mutation_filter_condition, 0, stop_value),
        depth = ifelse(clear_all_if_filter & mutation_filter_condition, 0, depth),
        A = ifelse(clear_all_if_filter & mutation_filter_condition, 0, A),
        C = ifelse(clear_all_if_filter & mutation_filter_condition, 0, C),
        G = ifelse(clear_all_if_filter & mutation_filter_condition, 0, G),
        `T` = ifelse(clear_all_if_filter & mutation_filter_condition, 0, `T`),
        N = ifelse(clear_all_if_filter & mutation_filter_condition, 0, N),
        P = ifelse(clear_all_if_filter & mutation_filter_condition, 0, P),
        PM = ifelse(clear_all_if_filter & mutation_filter_condition, 0, PM),
        S = ifelse(clear_all_if_filter & mutation_filter_condition, 0, S),
        SM = ifelse(clear_all_if_filter & mutation_filter_condition, 0, SM)
      ) %>%
      # Remove temporary column
      select(-mutation_filter_condition)
    cat("Mutation filtering applied. Threshold:", mutation_threshold, "Clear all data:", clear_all_if_filter, "\n")
  }
  
  ########################################################################
  # Apply stop filtering if requested
  ########################################################################
  if(filter_stop) {
    df <- df %>%
      mutate(
        # Store original stop condition to avoid circular dependency
        stop_filter_condition = normalized_stop_rate < stop_threshold,
        # Always clear main signal columns if normalized_stop_rate < stop_threshold
        mutation_rate = ifelse(stop_filter_condition, 0, mutation_rate),
        normalized_mutation_rate = ifelse(stop_filter_condition, 0, normalized_mutation_rate),
        stop_rate = ifelse(stop_filter_condition, 0, stop_rate),
        normalized_stop_rate = ifelse(stop_filter_condition, 0, normalized_stop_rate),
        signal_value = ifelse(stop_filter_condition, 0, signal_value),
        # Only clear technical columns if clear_all_if_filter = TRUE
        stop_value = ifelse(clear_all_if_filter & stop_filter_condition, 0, stop_value),
        depth = ifelse(clear_all_if_filter & stop_filter_condition, 0, depth),
        A = ifelse(clear_all_if_filter & stop_filter_condition, 0, A),
        C = ifelse(clear_all_if_filter & stop_filter_condition, 0, C),
        G = ifelse(clear_all_if_filter & stop_filter_condition, 0, G),
        `T` = ifelse(clear_all_if_filter & stop_filter_condition, 0, `T`),
        N = ifelse(clear_all_if_filter & stop_filter_condition, 0, N),
        P = ifelse(clear_all_if_filter & stop_filter_condition, 0, P),
        PM = ifelse(clear_all_if_filter & stop_filter_condition, 0, PM),
        S = ifelse(clear_all_if_filter & stop_filter_condition, 0, S),
        SM = ifelse(clear_all_if_filter & stop_filter_condition, 0, SM)
      ) %>%
      # Remove temporary column
      select(-stop_filter_condition)
    cat("Stop filtering applied. Threshold:", stop_threshold, "Clear all data:", clear_all_if_filter, "\n")
  }

  ########################################################################
  # Apply signal filtering if requested
  ########################################################################
  if(filter_signal) {
    df <- df %>%
      mutate(
        # Store original signal condition to avoid circular dependency
        signal_filter_condition = signal_value < signal_threshold,
        # Always clear main signal columns if signal_value < signal_threshold
        mutation_rate = ifelse(signal_filter_condition, 0, mutation_rate),
        normalized_mutation_rate = ifelse(signal_filter_condition, 0, normalized_mutation_rate),
        stop_rate = ifelse(signal_filter_condition, 0, stop_rate),
        normalized_stop_rate = ifelse(signal_filter_condition, 0, normalized_stop_rate),
        signal_value = ifelse(signal_filter_condition, 0, signal_value),
        # Only clear technical columns if clear_all_if_filter = TRUE
        stop_value = ifelse(clear_all_if_filter & signal_filter_condition, 0, stop_value),
        depth = ifelse(clear_all_if_filter & signal_filter_condition, 0, depth),
        A = ifelse(clear_all_if_filter & signal_filter_condition, 0, A),
        C = ifelse(clear_all_if_filter & signal_filter_condition, 0, C),
        G = ifelse(clear_all_if_filter & signal_filter_condition, 0, G),
        `T` = ifelse(clear_all_if_filter & signal_filter_condition, 0, `T`),
        N = ifelse(clear_all_if_filter & signal_filter_condition, 0, N),
        P = ifelse(clear_all_if_filter & signal_filter_condition, 0, P),
        PM = ifelse(clear_all_if_filter & signal_filter_condition, 0, PM),
        S = ifelse(clear_all_if_filter & signal_filter_condition, 0, S),
        SM = ifelse(clear_all_if_filter & signal_filter_condition, 0, SM)
      ) %>%
      # Remove temporary column
      select(-signal_filter_condition)
    cat("Signal filtering applied. Threshold:", signal_threshold, "Clear all data:", clear_all_if_filter, "\n")
  }
  
  ########################################################################
  # Apply final DHU site filtering if requested
  ########################################################################
  if(filter_dhu_site) {
    original_rows <- nrow(df)
    
    df <- df %>%
      dplyr::filter(ref_base == "T") %>%
      dplyr::filter(depth >= dhu_site_min_depth) %>%
      dplyr::filter(normalized_mutation_rate >= 0.05) %>%
      dplyr::filter(signal_value >= 0.08) %>%
      dplyr::filter(position > 10) %>%
      dplyr::filter(signal_value < 1)
    
    # Apply chromosome filtering if pattern is provided
    if (!is.null(chr_filter_pattern) && chr_filter_pattern != "") {
      df <- df %>%
        dplyr::filter(grepl(chr_filter_pattern, chr))
    }
    
    df <- df %>%
      dplyr::mutate(site_name = paste0(chr, "-", position))
    
    final_rows <- nrow(df)
    cat("DHU site filtering applied (min depth:", dhu_site_min_depth, "):", original_rows, "->", final_rows, "rows\n")
  }
  
  # Save processed data if requested
  if(save_data && !is.null(data_output_file)) {
    cat("Saving processed data to:", data_output_file, "\n")
    rio::export(df, data_output_file)
  }
  
  cat("DHU data processing complete!\n")
  
  return(df)
}

#' Process DHU sequencing data with background correction (enhanced)
#'
#' Supports multiple background files (averaged), T-to-C dominance calculation,
#' consecutive-T signal refinement, and secondary normalization of corrected rates.
#'
#' @param main_file_path Path to the main data file
#' @param background_file_path Path(s) to background data file(s). Can be a single
#'   path (character string) or a vector of paths; when multiple files are provided
#'   their rates are averaged before subtraction.
#' @param min_mean_depth Minimum mean depth threshold for chromosome filtering (default: 10)
#' @param mutation_threshold Threshold for mutation signal detection (default: 0.05)
#' @param stop_threshold Threshold for stop signal detection (default: 0.05)
#' @param filter_only_t Whether to filter only T reference bases (default: FALSE)
#' @param change_low_cov_to_0 Whether to apply low coverage filtering (default: FALSE)
#' @param site_depth_to_0 Depth threshold for low coverage filtering (default: 50)
#' @param save_data Whether to save processed data (default: FALSE)
#' @param data_output_file Output file path for processed data (optional)
#' @param match_by_chr Whether to match background data by chromosome (default: TRUE)
#' @param target_chr Target chromosome to keep after processing (default: NULL = keep all)
#' @param filter_mutation Whether to filter positions with low corrected mutation rates (default: FALSE)
#' @param filter_stop Whether to filter positions with low corrected stop rates (default: FALSE)
#' @param filter_signal Whether to filter positions with low corrected signal values (default: FALSE)
#' @param filter_t2c Whether to filter positions with low T-to-C dominance (default: FALSE)
#' @param signal_threshold Threshold for signal filtering (default: 0.05)
#' @param t2c_threshold Minimum T-to-C dominance ratio C/(A+G+C) (default: 0.66)
#' @param consecutive_t_correction Whether to refine signal for adjacent T sites with real signal (default: TRUE)
#' @param consecutive_t_mut_threshold Corrected mutation rate a T must exceed to qualify for consecutive-T refinement (default: 0.05)
#' @param clear_all_if_filter Whether to clear all data including technical columns when filtering (default: FALSE)
#' @param filter_dhu_site Whether to apply final DHU site filtering (default: FALSE)
#' @param chr_filter_pattern Pattern to filter chromosome names (default: NULL)
#' @param max_chr_length Final-step filter: keep only rows whose chr_length is
#'   below this value (default: 200, suited to tRNA; drops long rRNA/snRNA).
#'   Set to NULL to disable.
#' @param dhu_site_min_depth Minimum depth required for DHU site calls (default: 50)
#' @return Processed data frame with background-corrected rates, t2c_dominance,
#'   normalized_corrected_* columns, and original data
process_dhu_with_background <- function(main_file_path,
                                        background_file_path,
                                        min_mean_depth = 10,
                                        mutation_threshold = 0.05,
                                        stop_threshold = 0.05,
                                        filter_only_t = FALSE,
                                        change_low_cov_to_0 = FALSE,
                                        site_depth_to_0 = 50,
                                        save_data = FALSE,
                                        data_output_file = NULL,
                                        match_by_chr = TRUE,
                                        target_chr = NULL,
                                        filter_mutation = FALSE,
                                        filter_stop = FALSE,
                                        filter_signal = FALSE,
                                        filter_t2c = FALSE,
                                        signal_threshold = 0.05,
                                        t2c_threshold = 0.66,
                                        consecutive_t_correction = TRUE,
                                        consecutive_t_mut_threshold = 0.05,
                                        clear_all_if_filter = FALSE,
                                        filter_dhu_site = FALSE,
                                        chr_filter_pattern = NULL,
                                        max_chr_length = 200,
                                        dhu_site_min_depth = 50) {

  # Load required libraries
  require(dplyr)
  require(tidyr)
  require(rio)

  cat("=== process_dhu_with_background (enhanced) ===\n")
  cat("Main file:", main_file_path, "\n")
  cat("Background file(s):", paste(background_file_path, collapse = ", "), "\n")
  if (!is.null(target_chr)) cat("Target chr:", target_chr, "\n")
  cat("T2C dominance filtering:", ifelse(filter_t2c, paste("ENABLED (threshold:", t2c_threshold, ")"), "DISABLED"), "\n")
  cat("Consecutive-T correction:", ifelse(consecutive_t_correction,
      paste("ENABLED (mut >", consecutive_t_mut_threshold, "& dominance >=", t2c_threshold, ")"), "DISABLED"), "\n\n")

  ########################################################################
  # Process main file
  ########################################################################
  main_data <- process_dhu_data(
    file_path = main_file_path,
    min_mean_depth = min_mean_depth,
    mutation_threshold = mutation_threshold,
    stop_threshold = stop_threshold,
    filter_only_t = filter_only_t,
    change_low_cov_to_0 = change_low_cov_to_0,
    site_depth_to_0 = site_depth_to_0,
    save_data = FALSE,
    data_output_file = NULL,
    filter_mutation = FALSE,
    filter_stop = FALSE,
    filter_signal = FALSE,
    signal_threshold = signal_threshold,
    clear_all_if_filter = clear_all_if_filter
  )
  
  ########################################################################
  # Process background file(s) - support multiple, averaged
  ########################################################################
  bg_data_list <- lapply(background_file_path, function(bf) {
    bg <- process_dhu_data(
      file_path = bf,
      min_mean_depth = min_mean_depth,
      mutation_threshold = mutation_threshold,
      stop_threshold = stop_threshold,
      filter_only_t = filter_only_t,
      change_low_cov_to_0 = change_low_cov_to_0,
      site_depth_to_0 = site_depth_to_0,
      save_data = FALSE,
      data_output_file = NULL,
      filter_mutation = FALSE,
      filter_stop = FALSE,
      filter_signal = FALSE,
      signal_threshold = signal_threshold,
      clear_all_if_filter = clear_all_if_filter
    )
    
    # Select columns for background joining
    if (match_by_chr) {
      bg %>% select(chr, position,
                    bg_mutation_rate = mutation_rate,
                    bg_normalized_mutation_rate = normalized_mutation_rate,
                    bg_stop_rate = stop_rate,
                    bg_normalized_stop_rate = normalized_stop_rate,
                    bg_signal_value = signal_value)
    } else {
      bg %>% select(position,
                    bg_mutation_rate = mutation_rate,
                    bg_normalized_mutation_rate = normalized_mutation_rate,
                    bg_stop_rate = stop_rate,
                    bg_normalized_stop_rate = normalized_stop_rate,
                    bg_signal_value = signal_value)
    }
  })
  
  # Average background if multiple files
  if (match_by_chr) {
    join_keys <- c("chr", "position")
  } else {
    join_keys <- "position"
  }
  
  background_for_join <- dhu_prepare_background_for_join(bg_data_list, join_keys)
  if (length(bg_data_list) > 1) {
    cat("Averaged background from", length(bg_data_list), "files\n")
  }
  
  ########################################################################
  # Optional: filter to target_chr
  ########################################################################
  if (!is.null(target_chr)) {
    main_data <- main_data %>% filter(chr == target_chr)
    if ("chr" %in% names(background_for_join)) {
      background_for_join <- background_for_join %>% filter(chr == target_chr)
    }
    cat("Filtered to target_chr:", target_chr, "- main:", nrow(main_data), "rows; background:", nrow(background_for_join), "rows\n")
  }
  
  ########################################################################
  # Perform background correction
  ########################################################################
  corrected_data <- main_data %>%
    left_join(background_for_join, by = join_keys) %>%
    mutate(
      # Replace NA background values with 0
      bg_mutation_rate = ifelse(is.na(bg_mutation_rate), 0, bg_mutation_rate),
      bg_normalized_mutation_rate = ifelse(is.na(bg_normalized_mutation_rate), 0, bg_normalized_mutation_rate),
      bg_stop_rate = ifelse(is.na(bg_stop_rate), 0, bg_stop_rate),
      bg_normalized_stop_rate = ifelse(is.na(bg_normalized_stop_rate), 0, bg_normalized_stop_rate),
      bg_signal_value = ifelse(is.na(bg_signal_value), 0, bg_signal_value),
      
      # Calculate corrected rates (ensure non-negative values)
      correct_mutation = pmax(mutation_rate - bg_mutation_rate, 0),
      correct_normalized_mutation = pmax(normalized_mutation_rate - bg_normalized_mutation_rate, 0),
      correct_stop = pmax(stop_rate - bg_stop_rate, 0),
      correct_normalized_stop = pmax(normalized_stop_rate - bg_normalized_stop_rate, 0),
      
      # T-to-C dominance: C / (A + G + C)
      t2c_dominance = ifelse((A + G + C) > 0, C / (A + G + C), 0)
    )
  
  ########################################################################
  # Consecutive-T signal refinement
  ########################################################################
  if (consecutive_t_correction) {
    corrected_data <- corrected_data %>%
      arrange(chr, strand, position) %>%
      group_by(chr, strand) %>%
      mutate(
        cons_seq_len = max(position),
        cons_lead_S = lead(S, default = 0),
        cons_lead_SM = lead(SM, default = 0),
        # A T site qualifies if: ref==T, corrected mutation above threshold,
        # T-to-C dominated, has depth, and not in last 20 positions
        cons_self_qual = ref_base == "T" &
                         correct_mutation > consecutive_t_mut_threshold &
                         t2c_dominance >= t2c_threshold &
                         depth > 0 &
                         position <= (cons_seq_len - 20),
        cons_prev_qual = coalesce(lag(ref_base) == "T" &
                                  lag(correct_mutation) > consecutive_t_mut_threshold &
                                  lag(t2c_dominance) >= t2c_threshold, FALSE),
        cons_next_qual = coalesce(lead(ref_base) == "T" &
                                  lead(correct_mutation) > consecutive_t_mut_threshold &
                                  lead(t2c_dominance) >= t2c_threshold, FALSE),
        # Recalculate signal_value for consecutive-T pairs
        signal_value = case_when(
          # Latter T: subtract half of its own (S + SM)
          cons_self_qual & cons_prev_qual ~
            ((depth - (S + SM) / 2 - P) + cons_lead_S + cons_lead_SM) /
            (depth + cons_lead_S + cons_lead_SM),
          # Former T: count only half of the next site's (S + SM)
          cons_self_qual & cons_next_qual ~
            ((depth - P) + (cons_lead_S + cons_lead_SM) / 2) /
            (depth + (cons_lead_S + cons_lead_SM) / 2),
          # All other positions keep original signal_value
          TRUE ~ signal_value
        ),
        signal_value = pmax(signal_value, 0)
      ) %>%
      ungroup() %>%
      select(-cons_seq_len, -cons_lead_S, -cons_lead_SM,
             -cons_self_qual, -cons_prev_qual, -cons_next_qual)
    
    cat("Consecutive-T signal refinement applied\n")
  }
  
  ########################################################################
  # Background correction of signal (after possible consecutive-T refinement)
  ########################################################################
  corrected_data <- corrected_data %>%
    mutate(
      correct_signal = pmax(signal_value - bg_signal_value, 0)
    )
  
  ########################################################################
  # Secondary normalization of corrected rates
  ########################################################################
  corrected_data <- corrected_data %>%
    group_by(chr) %>%
    mutate(
      mean_correct_mutation = mean(correct_mutation, na.rm = TRUE),
      mean_correct_stop = mean(correct_stop, na.rm = TRUE),
      mean_correct_signal = mean(correct_signal, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      normalized_correct_mutation = pmax(correct_mutation - mean_correct_mutation, 0),
      normalized_correct_stop = pmax(correct_stop - mean_correct_stop, 0),
      normalized_correct_signal = pmax(correct_signal - mean_correct_signal, 0)
    )
  
  ########################################################################
  # Mark credible DHU sites (non-destructive; corrected values are preserved).
  # The is_dhu_site flag fully encodes the active threshold conditions, so no
  # signal columns are zeroed here. Set filter_dhu_site = TRUE to additionally
  # hard-filter rows down to DHU sites only.
  ########################################################################
  original_rows <- nrow(corrected_data)
  corrected_data <- dhu_apply_dhu_site_flag(
    corrected_data,
    filter_mutation = filter_mutation,
    filter_stop = filter_stop,
    filter_signal = filter_signal,
    filter_t2c = filter_t2c,
    mutation_threshold = mutation_threshold,
    stop_threshold = stop_threshold,
    signal_threshold = signal_threshold,
    t2c_threshold = t2c_threshold,
    dhu_site_min_depth = dhu_site_min_depth,
    chr_filter_pattern = chr_filter_pattern,
    hard_filter = filter_dhu_site
  )
  cat("DHU sites:", sum(corrected_data$is_dhu_site, na.rm = TRUE), "/", original_rows, "positions\n")
  if (filter_dhu_site) {
    cat("Hard-filtered to DHU sites:", original_rows, "->", nrow(corrected_data), "rows\n")
  }

  ########################################################################
  # Final filter: keep only short sequences (e.g. tRNA), drop long rRNA/snRNA
  ########################################################################
  if (!is.null(max_chr_length) && "chr_length" %in% names(corrected_data)) {
    rows_before <- nrow(corrected_data)
    corrected_data <- corrected_data %>% dplyr::filter(chr_length < max_chr_length)
    cat("chr_length filter (<", max_chr_length, "):", rows_before, "->", nrow(corrected_data), "rows\n")
  }

  ########################################################################
  # Summary statistics
  ########################################################################
  n_corrected_positions <- sum(corrected_data$correct_mutation > 0 | 
                               corrected_data$correct_stop > 0 | 
                               corrected_data$correct_signal > 0, na.rm = TRUE)
  cat("Positions with corrected signal:", n_corrected_positions, "\n")

  # Save processed data if requested
  if(save_data && !is.null(data_output_file)) {
    cat("Saving background-corrected data to:", data_output_file, "\n")
    rio::export(corrected_data, data_output_file)
  }

  cat("=== process_dhu_with_background completed ===\n")
  return(corrected_data)
}


#' Batch DHU analysis with background correction (enhanced wrapper)
#'
#' Processes multiple sample files against shared background file(s), combining
#' results into a single table. Replaces the 11+13 pipeline in a single call.
#'
#' @param files_to_process Vector of file paths to process
#' @param background_file_path Path(s) to background file(s), supports vector for averaging
#' @param output_dir Output directory for results (default: "analysis_output")
#' @param output_prefix Prefix for output file names (default: "batch_dhu")
#' @param target_chr Target chromosome to keep (default: NULL = keep all)
#' @param clear_first_positions Number of positions to set to 0 at chr start in filtered output (default: 5)
#' @param save_data Whether to write output files (default: TRUE)
#' @param min_mean_depth Minimum mean depth threshold (default: 10)
#' @param mutation_threshold Threshold for mutation filtering (default: 0.05)
#' @param stop_threshold Threshold for stop filtering (default: 0.05)
#' @param signal_threshold Threshold for signal filtering (default: 0.05)
#' @param t2c_threshold Minimum T-to-C dominance ratio (default: 0.66)
#' @param filter_only_t Whether to filter only T reference bases (default: FALSE)
#' @param change_low_cov_to_0 Whether to apply low coverage filtering (default: FALSE)
#' @param site_depth_to_0 Depth threshold for low coverage filtering (default: 50)
#' @param match_by_chr Whether to match background by chromosome (default: TRUE)
#' @param filter_mutation Whether to apply mutation filtering (default: FALSE)
#' @param filter_stop Whether to apply stop filtering (default: FALSE)
#' @param filter_signal Whether to apply signal filtering (default: FALSE)
#' @param filter_t2c Whether to apply T2C dominance filtering (default: FALSE)
#' @param consecutive_t_correction Whether to apply consecutive-T signal refinement (default: TRUE)
#' @param consecutive_t_mut_threshold Corrected mutation rate threshold for consecutive-T (default: 0.05)
#' @param clear_all_if_filter Whether to clear all columns when filtering (default: FALSE)
#' @param filter_dhu_site Whether to apply final DHU site filtering (default: FALSE)
#' @param chr_filter_pattern Pattern for chromosome filtering in filter_dhu_site (default: NULL)
#' @param max_chr_length Final-step filter: keep only rows whose chr_length is
#'   below this value (default: 200, suited to tRNA; drops long rRNA/snRNA).
#'   Set to NULL to disable.
#' @param dhu_site_min_depth Minimum depth required for DHU site calls (default: 50)
#' @return List with: processed_data (combined table), summary_report, parameters
batch_dhu_with_background <- function(files_to_process,
                                      background_file_path,
                                      output_dir = "analysis_output",
                                      output_prefix = "batch_dhu",
                                      target_chr = NULL,
                                      clear_first_positions = 5,
                                      save_data = TRUE,
                                      save_individual = TRUE,
                                      strip_from_name = NULL,
                                      min_mean_depth = 10,
                                      mutation_threshold = 0.05,
                                      stop_threshold = 0.05,
                                      signal_threshold = 0.05,
                                      t2c_threshold = 0.66,
                                      filter_only_t = FALSE,
                                      change_low_cov_to_0 = FALSE,
                                      site_depth_to_0 = 50,
                                      match_by_chr = TRUE,
                                      filter_mutation = FALSE,
                                      filter_stop = FALSE,
                                      filter_signal = FALSE,
                                      filter_t2c = FALSE,
                                      consecutive_t_correction = TRUE,
                                      consecutive_t_mut_threshold = 0.05,
                                      clear_all_if_filter = FALSE,
                                      filter_dhu_site = FALSE,
                                      chr_filter_pattern = NULL,
                                      max_chr_length = 200,
                                      dhu_site_min_depth = 50) {

  require(dplyr)
  require(tidyr)
  require(rio)
  
  cat("=== Batch DHU Analysis (Enhanced) ===\n")
  cat("Processing", length(files_to_process), "files\n")
  cat("Background file(s):", paste(background_file_path, collapse = ", "), "\n")
  if (!is.null(target_chr)) cat("Target chr:", target_chr, "\n")
  
  # If target_chr is specified, append chr name as subdirectory
  if (!is.null(target_chr)) {
    output_dir <- file.path(output_dir, target_chr)
  }
  
  cat("Output directory:", output_dir, "\n\n")
  
  # Create output directory
  if (save_data && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  ########################################################################
  # Pre-compute background (ONCE, then reuse for all samples)
  ########################################################################
  cat("--- Pre-computing background data (once) ---\n")
  
  if (match_by_chr) {
    join_keys <- c("chr", "position")
  } else {
    join_keys <- "position"
  }
  
  bg_data_list <- lapply(background_file_path, function(bf) {
    cat("  Processing background:", basename(bf), "\n")
    bg <- process_dhu_data(
      file_path = bf,
      min_mean_depth = min_mean_depth,
      mutation_threshold = mutation_threshold,
      stop_threshold = stop_threshold,
      filter_only_t = filter_only_t,
      change_low_cov_to_0 = change_low_cov_to_0,
      site_depth_to_0 = site_depth_to_0,
      save_data = FALSE,
      data_output_file = NULL,
      filter_mutation = FALSE,
      filter_stop = FALSE,
      filter_signal = FALSE,
      signal_threshold = signal_threshold,
      clear_all_if_filter = FALSE
    )
    
    if (match_by_chr) {
      bg %>% select(chr, position,
                    bg_mutation_rate = mutation_rate,
                    bg_normalized_mutation_rate = normalized_mutation_rate,
                    bg_stop_rate = stop_rate,
                    bg_normalized_stop_rate = normalized_stop_rate,
                    bg_signal_value = signal_value)
    } else {
      bg %>% select(position,
                    bg_mutation_rate = mutation_rate,
                    bg_normalized_mutation_rate = normalized_mutation_rate,
                    bg_stop_rate = stop_rate,
                    bg_normalized_stop_rate = normalized_stop_rate,
                    bg_signal_value = signal_value)
    }
  })
  
  # Average background by join key. This also collapses duplicate positions
  # when match_by_chr = FALSE, even for a single background file.
  background_for_join <- dhu_prepare_background_for_join(bg_data_list, join_keys)
  if (length(bg_data_list) > 1) {
    cat("  Averaged background from", length(bg_data_list), "files\n")
  }
  
  # Optional: filter background to target_chr
  if (!is.null(target_chr) && "chr" %in% names(background_for_join)) {
    background_for_join <- background_for_join %>% filter(chr == target_chr)
  }
  
  cat("  Background ready:", nrow(background_for_join), "positions\n\n")
  
  ########################################################################
  # Process each sample file (reusing pre-computed background)
  ########################################################################
  processed_results <- list()
  processed_sample_names <- character()
  
  for (i in seq_along(files_to_process)) {
    file_path <- files_to_process[i]
    sample_name <- gsub("-merged\\.txt$", "", basename(file_path))
    # Strip unwanted patterns from sample name
    if (!is.null(strip_from_name)) {
      for (pat in strip_from_name) {
        sample_name <- gsub(pat, "", sample_name, fixed = TRUE)
      }
    }
    cat("--- Processing file", i, "of", length(files_to_process), ":", basename(file_path), "---\n")
    
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path, ". Skipping...")
      next
    }
    processed_sample_names <- c(processed_sample_names, sample_name)
    
    # Process main file
    main_data <- process_dhu_data(
      file_path = file_path,
      min_mean_depth = min_mean_depth,
      mutation_threshold = mutation_threshold,
      stop_threshold = stop_threshold,
      filter_only_t = filter_only_t,
      change_low_cov_to_0 = change_low_cov_to_0,
      site_depth_to_0 = site_depth_to_0,
      save_data = FALSE,
      data_output_file = NULL,
      filter_mutation = FALSE,
      filter_stop = FALSE,
      filter_signal = FALSE,
      signal_threshold = signal_threshold,
      clear_all_if_filter = FALSE
    )
    
    # Filter to target_chr if specified
    if (!is.null(target_chr)) {
      main_data <- main_data %>% filter(chr == target_chr)
    }
    
    # Background correction
    corrected_data <- main_data %>%
      left_join(background_for_join, by = join_keys) %>%
      mutate(
        bg_mutation_rate = ifelse(is.na(bg_mutation_rate), 0, bg_mutation_rate),
        bg_normalized_mutation_rate = ifelse(is.na(bg_normalized_mutation_rate), 0, bg_normalized_mutation_rate),
        bg_stop_rate = ifelse(is.na(bg_stop_rate), 0, bg_stop_rate),
        bg_normalized_stop_rate = ifelse(is.na(bg_normalized_stop_rate), 0, bg_normalized_stop_rate),
        bg_signal_value = ifelse(is.na(bg_signal_value), 0, bg_signal_value),
        correct_mutation = pmax(mutation_rate - bg_mutation_rate, 0),
        correct_normalized_mutation = pmax(normalized_mutation_rate - bg_normalized_mutation_rate, 0),
        correct_stop = pmax(stop_rate - bg_stop_rate, 0),
        correct_normalized_stop = pmax(normalized_stop_rate - bg_normalized_stop_rate, 0),
        t2c_dominance = ifelse((A + G + C) > 0, C / (A + G + C), 0)
      )
    
    # Consecutive-T signal refinement
    if (consecutive_t_correction) {
      corrected_data <- corrected_data %>%
        arrange(chr, strand, position) %>%
        group_by(chr, strand) %>%
        mutate(
          cons_seq_len = max(position),
          cons_lead_S = lead(S, default = 0),
          cons_lead_SM = lead(SM, default = 0),
          cons_self_qual = ref_base == "T" &
                           correct_mutation > consecutive_t_mut_threshold &
                           t2c_dominance >= t2c_threshold &
                           depth > 0 &
                           position <= (cons_seq_len - 20),
          cons_prev_qual = coalesce(lag(ref_base) == "T" &
                                    lag(correct_mutation) > consecutive_t_mut_threshold &
                                    lag(t2c_dominance) >= t2c_threshold, FALSE),
          cons_next_qual = coalesce(lead(ref_base) == "T" &
                                    lead(correct_mutation) > consecutive_t_mut_threshold &
                                    lead(t2c_dominance) >= t2c_threshold, FALSE),
          signal_value = case_when(
            cons_self_qual & cons_prev_qual ~
              ((depth - (S + SM) / 2 - P) + cons_lead_S + cons_lead_SM) /
              (depth + cons_lead_S + cons_lead_SM),
            cons_self_qual & cons_next_qual ~
              ((depth - P) + (cons_lead_S + cons_lead_SM) / 2) /
              (depth + (cons_lead_S + cons_lead_SM) / 2),
            TRUE ~ signal_value
          ),
          signal_value = pmax(signal_value, 0)
        ) %>%
        ungroup() %>%
        select(-cons_seq_len, -cons_lead_S, -cons_lead_SM,
               -cons_self_qual, -cons_prev_qual, -cons_next_qual)
    }
    
    # Background correction of signal (after consecutive-T refinement)
    corrected_data <- corrected_data %>%
      mutate(correct_signal = pmax(signal_value - bg_signal_value, 0))
    
    # Secondary normalization
    corrected_data <- corrected_data %>%
      group_by(chr) %>%
      mutate(
        mean_correct_mutation = mean(correct_mutation, na.rm = TRUE),
        mean_correct_stop = mean(correct_stop, na.rm = TRUE),
        mean_correct_signal = mean(correct_signal, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      mutate(
        normalized_correct_mutation = pmax(correct_mutation - mean_correct_mutation, 0),
        normalized_correct_stop = pmax(correct_stop - mean_correct_stop, 0),
        normalized_correct_signal = pmax(correct_signal - mean_correct_signal, 0)
      )
    
    # Mark credible DHU sites non-destructively (corrected values preserved).
    # is_dhu_site fully encodes the active thresholds; set filter_dhu_site = TRUE
    # to additionally hard-filter the table down to DHU sites only.
    corrected_data <- dhu_apply_dhu_site_flag(
      corrected_data,
      filter_mutation = filter_mutation,
      filter_stop = filter_stop,
      filter_signal = filter_signal,
      filter_t2c = filter_t2c,
      mutation_threshold = mutation_threshold,
      stop_threshold = stop_threshold,
      signal_threshold = signal_threshold,
      t2c_threshold = t2c_threshold,
      dhu_site_min_depth = dhu_site_min_depth,
      chr_filter_pattern = chr_filter_pattern,
      hard_filter = filter_dhu_site
    )

    # Final filter: keep only short sequences (e.g. tRNA), drop long rRNA/snRNA
    if (!is.null(max_chr_length) && "chr_length" %in% names(corrected_data)) {
      corrected_data <- corrected_data %>% dplyr::filter(chr_length < max_chr_length)
    }

    cat("  DHU sites:", sum(corrected_data$is_dhu_site, na.rm = TRUE), "/", nrow(corrected_data), "positions\n")

    # Add sample name
    corrected_data$sample_name <- sample_name
    
    # Save individual file result
    if (save_data && save_individual) {
      individual_file <- file.path(output_dir, paste0(sample_name, ".csv"))
      write.csv(corrected_data, file = individual_file, row.names = FALSE)
      cat("  Saved individual:", basename(individual_file), "\n")
    }
    
    processed_results[[i]] <- corrected_data
    cat("  Done:", nrow(corrected_data), "rows\n")
  }
  
  # Combine all results
  processed_results <- Filter(Negate(is.null), processed_results)
  if (length(processed_results) == 0) {
    stop("No input files were processed successfully.")
  }
  
  all_results <- bind_rows(processed_results)
  cat("\n=== Combined", length(processed_results), "samples,", nrow(all_results), "total rows ===\n")
  
  # Build summary report (based on is_dhu_site flag)
  summary_source <- dhu_zero_first_positions(all_results, clear_first_positions)
  summary_report <- summary_source %>%
    group_by(sample_name) %>%
    summarise(
      total_positions = n(),
      dhu_sites = sum(is_dhu_site, na.rm = TRUE),
      max_correct_mutation = dhu_safe_max(correct_mutation[is_dhu_site]),
      max_correct_stop = dhu_safe_max(correct_stop[is_dhu_site]),
      max_correct_signal = dhu_safe_max(correct_signal[is_dhu_site]),
      max_normalized_correct_mutation = dhu_safe_max(normalized_correct_mutation[is_dhu_site]),
      max_normalized_correct_stop = dhu_safe_max(normalized_correct_stop[is_dhu_site]),
      max_normalized_correct_signal = dhu_safe_max(normalized_correct_signal[is_dhu_site]),
      max_t2c_dominance = dhu_safe_max(t2c_dominance[is_dhu_site]),
      .groups = 'drop'
    )
  
  summary_report <- data.frame(sample_name = unique(processed_sample_names)) %>%
    dplyr::left_join(summary_report, by = "sample_name") %>%
    dplyr::mutate(
      total_positions = ifelse(is.na(total_positions), 0L, total_positions),
      dhu_sites = ifelse(is.na(dhu_sites), 0L, dhu_sites)
    )
  
  # Write combined outputs
  if (save_data) {
    raw_csv <- file.path(output_dir, paste0(output_prefix, "_all_processed_data.csv"))
    write.csv(all_results, file = raw_csv, row.names = FALSE)
    write.table(all_results,
                file = file.path(output_dir, paste0(output_prefix, "_all_processed_data.txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    filtered_data <- dhu_zero_first_positions(all_results, clear_first_positions)
    filtered_csv <- file.path(
      output_dir,
      paste0(output_prefix, "_filtered_data_pos", clear_first_positions + 1, "plus.csv")
    )
    write.csv(filtered_data, file = filtered_csv, row.names = FALSE)

    # DHU sites only
    dhu_only <- all_results %>% dplyr::filter(is_dhu_site)
    write.csv(dhu_only,
              file = file.path(output_dir, paste0(output_prefix, "_dhu_sites.csv")),
              row.names = FALSE)

    write.csv(summary_report,
              file = file.path(output_dir, paste0(output_prefix, "_summary_report.csv")),
              row.names = FALSE)

    cat("\nOutput files written to:", output_dir, "\n")
    cat("  Full table:", basename(raw_csv), "(", nrow(all_results), "rows )\n")
    cat("  Filtered table:", basename(filtered_csv), "(", nrow(filtered_data), "rows )\n")
    cat("  DHU sites:", nrow(dhu_only), "rows\n")
    cat("  Individual files:", length(processed_results), "x {sample_name}.csv\n")
  }
  
  cat("=== Batch DHU Analysis Completed ===\n")
  
  return(list(
    processed_data = all_results,
    dhu_sites = all_results %>% dplyr::filter(is_dhu_site),
    summary_report = summary_report,
    parameters = list(
      files_processed = files_to_process,
      background_file_path = background_file_path,
      target_chr = target_chr,
      output_prefix = output_prefix,
      output_dir = output_dir,
      clear_first_positions = clear_first_positions,
      mutation_threshold = mutation_threshold,
      stop_threshold = stop_threshold,
      signal_threshold = signal_threshold,
      t2c_threshold = t2c_threshold,
      filter_mutation = filter_mutation,
      filter_stop = filter_stop,
      filter_signal = filter_signal,
      filter_t2c = filter_t2c,
      filter_dhu_site = filter_dhu_site,
      dhu_site_min_depth = dhu_site_min_depth,
      max_chr_length = max_chr_length,
      consecutive_t_correction = consecutive_t_correction,
      consecutive_t_mut_threshold = consecutive_t_mut_threshold
    )
  ))
}
