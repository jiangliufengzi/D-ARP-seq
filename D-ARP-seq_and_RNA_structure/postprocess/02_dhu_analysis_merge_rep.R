# DHU 分析重复样本合并 - 步骤 2
# 作者：为 DHU 分析流程生成
# 
# 包含合并两个经 dhu_analysis_step1.R 处理的重复样本文件的函数
# 主函数：merge_dhu_replicates() - 合并两个经背景校正的 DHU 重复样本

#' 合并两个经背景校正的 DHU 重复样本文件
#'
#' 该函数合并两个已通过 01_dhu_analysis.R 处理并完成背景校正的重复样本文件。
#' 对于两个重复样本均存在的位点，数值信号列取平均值；结构列(chr_length 等)
#' 取 rep1 值；is_dhu_site 按 dhu_site_rule 合并，site_name 据合并后的
#' is_dhu_site 重新生成以保持一致，sample_name 设为单一合并名。
#' 对于仅出现在一个重复样本中的位点，行为由 keep_in_one_file 参数控制。
#'
#' @param replicate1_file 第一个重复样本文件路径（已背景校正）
#' @param replicate2_file 第二个重复样本文件路径（已背景校正）
#' @param keep_in_one_file 逻辑值。若为 TRUE，保留仅出现在一个重复样本中的位点并使用原始值；若为 FALSE，则移除此类位点（默认：TRUE）
#' @param output_file 可选，合并后数据的保存路径（默认：NULL，不保存）
#' @param replicate1_name 重复样本 1 的名称标识（默认："rep1"）
#' @param replicate2_name 重复样本 2 的名称标识（默认："rep2"）
#' @param dhu_site_rule is_dhu_site 的合并规则："and"=两重复都为TRUE才算(默认,
#'   最严谨/可重复)，"or"=任一为TRUE即算，"rep1"=沿用 rep1 的判定。单重复位点
#'   在 "and" 下恒为 FALSE(无法判定可重复性)
#' @param merged_sample_name 合并后 sample_name 的值（默认 NULL：自动从 rep1 的
#'   sample_name 去掉 rep1/rep2 标记推断）
#' @param take_first_cols 在两重复间应相同、直接取 rep1 而非平均的数值列
#'   （默认 c("chr_length","seq_length")）
#' @param verbose 逻辑值。若为 TRUE，打印进度信息（默认：TRUE）
#' @return 合并后的重复样本数据框
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
  
  # 加载所需 R 包
  require(dplyr)
  require(tidyr)
  require(rio)
  
  if(verbose) {
    cat("=== DHU Replicate Merger (Step 2) ===\n")
    cat("Replicate 1:", replicate1_file, "\n")
    cat("Replicate 2:", replicate2_file, "\n")
    cat("Keep single-replicate positions:", keep_in_one_file, "\n")
  }
  
  # 读取重复样本文件
  if(verbose) cat("Loading replicate files...\n")
  rep1_data <- normalize_replicate_import(rio::import(replicate1_file))
  rep2_data <- normalize_replicate_import(rio::import(replicate2_file))
  
  if(verbose) {
    cat("Replicate 1 data:", nrow(rep1_data), "rows,", ncol(rep1_data), "columns\n")
    cat("Replicate 2 data:", nrow(rep2_data), "rows,", ncol(rep2_data), "columns\n")
  }
  
  # 验证数据结构
  validate_replicate_data(rep1_data, rep2_data, verbose)
  
  # 添加重复样本来源标识
  rep1_data$replicate_source <- replicate1_name
  rep2_data$replicate_source <- replicate2_name
  
  # 定义用于连接的键列（chr + position）
  join_keys <- c("chr", "position")

  # 这几列在合并后单独处理（不走通用平均/rep1优先逻辑）
  special_cols <- intersect(c("is_dhu_site", "site_name", "sample_name"), names(rep1_data))

  # 识别需要取平均值的数值列（排除键列、来源列、结构列、特殊列）
  numeric_cols <- identify_numeric_columns(
    rep1_data,
    exclude_cols = c(join_keys, "replicate_source", take_first_cols, special_cols)
  )

  if(verbose) {
    cat("Joining by:", paste(join_keys, collapse = ", "), "\n")
    cat("Numeric columns to average:", length(numeric_cols), "columns\n")
    cat("is_dhu_site merge rule:", dhu_site_rule, "\n")
  }

  # 执行全连接以保留所有位点
  merged_data <- full_join(rep1_data, rep2_data,
                          by = join_keys,
                          suffix = c(".rep1", ".rep2"))

  # 计算数据来源信息
  merged_data <- merged_data %>%
    mutate(
      # 确定数据来源
      data_source = case_when(
        !is.na(replicate_source.rep1) & !is.na(replicate_source.rep2) ~ "both_replicates",
        !is.na(replicate_source.rep1) & is.na(replicate_source.rep2) ~ replicate1_name,
        is.na(replicate_source.rep1) & !is.na(replicate_source.rep2) ~ replicate2_name,
        TRUE ~ "unknown"
      )
    )

  # (1) 数值信号列：两重复都有则取平均，单重复按 keep_in_one_file 处理
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

  # (2) 结构列（两重复应相同）：直接取 rep1，缺失则取 rep2，不平均
  for(col in intersect(take_first_cols, names(rep1_data))) {
    col_rep1 <- paste0(col, ".rep1")
    col_rep2 <- paste0(col, ".rep2")
    if(col_rep1 %in% names(merged_data) && col_rep2 %in% names(merged_data)) {
      merged_data[[col]] <- coalesce(merged_data[[col_rep1]], merged_data[[col_rep2]])
      merged_data[[col_rep1]] <- NULL
      merged_data[[col_rep2]] <- NULL
    }
  }

  # (3) 其余元数据非数值列（ref_base/strand/sprinzl_position/*_yes_or_no）：rep1 优先
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

  # (4) is_dhu_site：按规则合并（单重复缺失视为 FALSE）
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

  # (5) site_name：据合并后的 is_dhu_site 重新生成，保证一致
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

  # (6) sample_name：合并为单一名称
  for(sn in c("sample_name.rep1", "sample_name.rep2")) {
    if(sn %in% names(merged_data)) merged_data[[sn]] <- NULL
  }
  if("sample_name" %in% names(rep1_data)) {
    if(is.null(merged_sample_name)) {
      s1 <- as.character(rep1_data$sample_name[1])
      # 去掉以分隔符开头的 rep1/rep2 标记及其后缀（要求 rep 前有 . _ - 分隔符，
      # 避免误删名字中间含 "rep" 的部分）
      merged_sample_name <- sub("[._-]rep[12].*$", "", s1)
      if(is.na(merged_sample_name) || merged_sample_name == "") merged_sample_name <- s1
    }
    merged_data$sample_name <- merged_sample_name
  }

  # 清理重复样本来源列
  merged_data$replicate_source.rep1 <- NULL
  merged_data$replicate_source.rep2 <- NULL
  
  # 根据 keep_in_one_file 参数过滤数据
  if(!keep_in_one_file) {
    original_rows <- nrow(merged_data)
    merged_data <- merged_data %>%
      filter(data_source == "both_replicates")
    if(verbose) {
      cat("Removed", original_rows - nrow(merged_data), "positions present in only one replicate\n")
    }
  }

  # 列顺序：chr, position, sprinzl_position 居前三列，其余按 rep1 原顺序排列
  if("sprinzl_position" %in% names(merged_data)) {
    lead_cols <- intersect(c("chr", "position", "sprinzl_position"), names(merged_data))
    tail_cols <- setdiff(names(merged_data), lead_cols)
    ordered_tail <- c(
      intersect(names(rep1_data), tail_cols),
      setdiff(tail_cols, names(rep1_data))
    )
    merged_data <- merged_data[, c(lead_cols, ordered_tail), drop = FALSE]
  }
  
  # 添加汇总统计
  summary_stats <- merged_data %>%
    group_by(data_source) %>%
    summarise(count = n(), .groups = 'drop')
  
  if(verbose) {
    cat("\nMerge summary:\n")
    print(summary_stats)
    cat("Final merged data:", nrow(merged_data), "rows,", ncol(merged_data), "columns\n")
  }
  
  # 按需保存数据
  if(!is.null(output_file)) {
    if(verbose) cat("Saving merged data to:", output_file, "\n")
    rio::export(merged_data, output_file)
  }
  
  if(verbose) cat("DHU replicate merging completed successfully!\n")
  
  return(merged_data)
}

#' 验证重复样本数据结构及兼容性
#'
#' @param rep1_data 第一个重复样本数据框
#' @param rep2_data 第二个重复样本数据框
#' @param verbose 是否输出进度信息
validate_replicate_data <- function(rep1_data, rep2_data, verbose = TRUE) {
  
  # 检查必需列是否存在
  required_cols <- c("chr", "position")
  
  missing_rep1 <- setdiff(required_cols, names(rep1_data))
  missing_rep2 <- setdiff(required_cols, names(rep2_data))
  
  if(length(missing_rep1) > 0) {
    stop("Replicate 1 missing required columns: ", paste(missing_rep1, collapse = ", "))
  }
  
  if(length(missing_rep2) > 0) {
    stop("Replicate 2 missing required columns: ", paste(missing_rep2, collapse = ", "))
  }
  
  # 检查列兼容性
  common_cols <- intersect(names(rep1_data), names(rep2_data))
  if(verbose) {
    cat("Common columns between replicates:", length(common_cols), "\n")
  }
  
  # 对列差异发出警告
  rep1_only <- setdiff(names(rep1_data), names(rep2_data))
  rep2_only <- setdiff(names(rep2_data), names(rep1_data))
  
  if(length(rep1_only) > 0 && verbose) {
    cat("Columns only in replicate 1:", paste(rep1_only, collapse = ", "), "\n")
  }
  
  if(length(rep2_only) > 0 && verbose) {
    cat("Columns only in replicate 2:", paste(rep2_only, collapse = ", "), "\n")
  }
}

#' 识别需要取平均值的数值列
#'
#' @param data 待分析的数据框
#' @param exclude_cols 不参与取平均值的列
#' @return 数值列名称向量
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
  
  # 获取所有数值列
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  
  # 排除指定列
  numeric_cols <- setdiff(numeric_cols, exclude_cols)
  
  return(numeric_cols)
}

#' 计算一对重复样本列的合并值
#'
#' @param val1 重复样本 1 的值
#' @param val2 重复样本 2 的值
#' @param keep_single 是否保留仅出现在一个重复样本中的值
#' @return 合并后的值向量
calculate_merged_value <- function(val1, val2, keep_single = TRUE) {
  
  # 初始化结果向量
  result <- rep(NA, length(val1))
  
  # 两个值均存在：计算平均值
  both_present <- !is.na(val1) & !is.na(val2)
  result[both_present] <- (val1[both_present] + val2[both_present]) / 2
  
  if(keep_single) {
    # 仅 val1 存在
    only_val1 <- !is.na(val1) & is.na(val2)
    result[only_val1] <- val1[only_val1]
    
    # 仅 val2 存在
    only_val2 <- is.na(val1) & !is.na(val2)
    result[only_val2] <- val2[only_val2]
  }
  
  return(result)
}

#' 批量合并多对重复样本
#'
#' @param replicate_pairs 列表的列表，每项包含 replicate1_file 和 replicate2_file 路径
#' @param output_dir 保存合并文件的目录（可选）
#' @param keep_in_one_file 传递给 merge_dhu_replicates 的参数
#' @param verbose 是否输出进度信息
#' @return 合并后数据框的列表
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
    
    # 若指定了输出目录，则生成输出文件名
    output_file <- NULL
    if(!is.null(output_dir)) {
      base_name <- paste0("merged_pair_", i, "_", Sys.Date(), ".csv")
      output_file <- file.path(output_dir, base_name)
    }
    
    # 合并当前这对重复样本
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
