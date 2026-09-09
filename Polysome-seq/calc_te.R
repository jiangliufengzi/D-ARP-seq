#!/usr/bin/env Rscript
# Manual translational efficiency. Matches li_pipeline 02_calc_TE.py /
# 02_qc_te/06_normalization.R apply_to_te path:
#   classic TPM, then translation-side TPM x te_ratio(Q), RNA/Input Q = 1
#   TE = (TPM_translation * Q + pseudo) / (TPM_abundance + pseudo)
#   both sides < min_tpm (after Q) -> NA
# Spike-in table must contain sample + te_ratio (not library Y/geomean(Y)).
#
# Usage:
#   Rscript calc_te.R --mode polysome --counts count_matrix.txt --output-dir 07.TE \
#       [--spikein-table spikein_size_factors.tsv] [--sample-sheet samples.tsv]
#   Rscript calc_te.R --mode ribo --ribo-counts ribo.txt --rna-counts rna.txt \
#       --lengths lengths.tsv --sample-sheet samples.tsv --output-dir 07.TE \
#       [--spikein-table spikein_size_factors.tsv]

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL, required = FALSE) {
  hit <- which(args == flag)
  if (!length(hit)) {
    if (required) stop("missing required argument: ", flag)
    return(default)
  }
  if (hit[1] == length(args)) stop("argument ", flag, " needs a value")
  args[hit[1] + 1]
}

clean_sample <- function(x) {
  x <- basename(as.character(x))
  x <- sub("\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("_Aligned\\.sortedByCoord\\.out$", "", x)
  x <- sub("_Aligned\\.toTranscriptome\\.out$", "", x)
  x <- sub("[-_]merged-filtered$", "", x)
  x <- sub("[-_]merged$", "", x)
  x <- sub("[-_]filtered$", "", x)
  x <- sub("[-_]unmapped$", "", x)
  x
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("file not found: ", path)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "#")
}

lookup_sf <- function(sf, sample) {
  if (is.null(sf)) return(NA_real_)
  key <- clean_sample(sample)
  if (key %in% names(sf)) return(unname(sf[[key]]))
  hit <- which(vapply(names(sf), function(k) grepl(k, key, fixed = TRUE), logical(1)))
  if (length(hit) == 1) return(unname(sf[[hit]]))
  hit <- which(vapply(names(sf), function(k) grepl(key, k, fixed = TRUE), logical(1)))
  if (length(hit) == 1) return(unname(sf[[hit]]))
  NA_real_
}

load_te_ratio <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) {
    warning("spike-in table not found, TE uses Q=1: ", path)
    return(NULL)
  }
  df <- read_tsv(path)
  names(df) <- tolower(names(df))
  if (!"sample" %in% names(df)) stop("spike-in table needs a sample column: ", path)
  ratio_col <- intersect(c("te_ratio", "ratio_poly_hy", "ratio"), names(df))
  if (!length(ratio_col)) {
    message(
      "[WARN] spike-in table has no te_ratio; manual TE will use Q=1. ",
      "Library size_factor=Y/geomean(Y) is not the original Q."
    )
    return(NULL)
  }
  tr <- as.numeric(df[[ratio_col[1]]])
  names(tr) <- clean_sample(df$sample)
  tr <- tr[is.finite(tr) & tr > 0]
  if (!length(tr)) stop("no positive te_ratio values in: ", path)
  message("[INFO] loaded ", length(tr), " te_ratio (Q) values from ", path, " column ", ratio_col[1])
  tr
}

tpm_from_counts <- function(counts, lengths) {
  id <- counts[[1]]
  mat <- as.matrix(counts[-1])
  storage.mode(mat) <- "double"
  len <- as.numeric(lengths[match(id, names(lengths))])
  len_kb <- len / 1000
  len_kb[!is.finite(len_kb) | len_kb <= 0] <- NA_real_
  rpk <- mat / len_kb
  rpk[!is.finite(rpk)] <- 0
  scale <- colSums(rpk, na.rm = TRUE) / 1e6
  scale[!is.finite(scale) | scale <= 0] <- 1
  tpm <- sweep(rpk, 2, scale, "/")
  data.frame(id = id, tpm, check.names = FALSE, stringsAsFactors = FALSE)
}

bh_adjust <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (!any(ok)) return(out)
  pv <- p[ok]
  n <- length(pv)
  o <- order(pv)
  adj <- pv[o] * n / seq_len(n)
  adj <- pmin(1, rev(cummin(rev(adj))))
  out_ok <- numeric(n)
  out_ok[o] <- adj
  out[ok] <- out_ok
  out
}

welch_rows <- function(ctrl, trt, min_n = 2) {
  n <- nrow(ctrl)
  mean_c <- mean_t <- tstat <- pval <- rep(NA_real_, n)
  n_c <- n_t <- integer(n)
  for (i in seq_len(n)) {
    x <- as.numeric(ctrl[i, ])
    y <- as.numeric(trt[i, ])
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    n_c[i] <- length(x)
    n_t[i] <- length(y)
    if (length(x)) mean_c[i] <- mean(x)
    if (length(y)) mean_t[i] <- mean(y)
    if (length(x) < min_n || length(y) < min_n) next
    if (stats::sd(c(x, y)) == 0) next
    res <- tryCatch(stats::t.test(y, x, var.equal = FALSE), error = function(e) NULL)
    if (!is.null(res)) {
      tstat[i] <- unname(res$statistic)
      pval[i] <- res$p.value
    }
  }
  data.frame(
    n_control = n_c,
    n_treatment = n_t,
    mean_log2TE_control = mean_c,
    mean_log2TE_treatment = mean_t,
    log2FoldChange = mean_t - mean_c,
    t_stat = tstat,
    pvalue = pval,
    padj = bh_adjust(pval),
    stringsAsFactors = FALSE
  )
}

write_te <- function(id, pairs, tpm_tr, tpm_ab, q, pseudo, min_tpm, out_dir, id_name) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  te <- data.frame(id, stringsAsFactors = FALSE, check.names = FALSE)
  names(te)[1] <- id_name
  summary_rows <- list()
  for (i in seq_len(nrow(pairs))) {
    pr <- pairs[i, ]
    tr <- as.numeric(tpm_tr[[pr$translation]])
    ab <- as.numeric(tpm_ab[[pr$abundance]])
    tr_corr <- tr * q[i]
    both_low <- (tr_corr < min_tpm) & (ab < min_tpm)
    ratio <- (tr_corr + pseudo) / (ab + pseudo)
    log2te <- log2(ratio)
    ratio[both_low] <- NA_real_
    log2te[both_low] <- NA_real_
    if (identical(id_name, "gene")) {
      te[[paste0("TPM_input_", pr$pair)]] <- ab
      te[[paste0("TPM_poly_", pr$pair)]] <- tr_corr
      te[[paste0("TE_", pr$pair)]] <- ratio
      te[[paste0("log2TE_", pr$pair)]] <- log2te
    } else {
      te[[paste0(pr$pair, "_log2TE")]] <- log2te
    }
    vals <- log2te[is.finite(log2te)]
    summary_rows[[i]] <- data.frame(
      pair = pr$pair,
      group = pr$group,
      translation_sample = pr$translation,
      abundance_sample = pr$abundance,
      Q = q[i],
      n_genes = length(vals),
      n_filtered = sum(both_low),
      median_log2TE = if (length(vals)) stats::median(vals) else NA_real_,
      mean_log2TE = if (length(vals)) mean(vals) else NA_real_,
      sd_log2TE = if (length(vals) > 1) stats::sd(vals) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  utils::write.table(te, file.path(out_dir, "TE_matrix.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  log2_cols <- grep("log2TE", names(te), value = TRUE)
  log2_only <- te[, c(id_name, log2_cols), drop = FALSE]
  utils::write.table(log2_only, file.path(out_dir, "te_log2_per_pair.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  summary <- do.call(rbind, summary_rows)
  utils::write.table(summary, file.path(out_dir, "TE_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(pairs, file.path(out_dir, "te_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  groups <- unique(pairs$group)
  if (length(groups) == 2) {
    g1 <- groups[1]
    g2 <- groups[2]
    c1 <- intersect(
      c(paste0("log2TE_", pairs$pair[pairs$group == g1]),
        paste0(pairs$pair[pairs$group == g1], "_log2TE")),
      names(te)
    )
    c2 <- intersect(
      c(paste0("log2TE_", pairs$pair[pairs$group == g2]),
        paste0(pairs$pair[pairs$group == g2], "_log2TE")),
      names(te)
    )
    if (length(c1) >= 2 && length(c2) >= 2) {
      stats <- welch_rows(as.matrix(te[, c1, drop = FALSE]), as.matrix(te[, c2, drop = FALSE]))
      stats <- cbind(setNames(data.frame(te[[1]], stringsAsFactors = FALSE), id_name), stats)
      names(stats)[names(stats) == "mean_log2TE_control"] <- paste0("mean_log2TE_", g1)
      names(stats)[names(stats) == "mean_log2TE_treatment"] <- paste0("mean_log2TE_", g2)
      stats$control_group <- g1
      stats$treatment_group <- g2
      utils::write.table(stats, file.path(out_dir, "TE_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  invisible(te)
}

parse_featurecounts <- function(path) {
  df <- read_tsv(path)
  names(df)[names(df) == "Geneid"] <- "id"
  if (!"id" %in% names(df)) names(df)[1] <- "id"
  if (!"Length" %in% names(df)) stop("featureCounts matrix needs a Length column: ", path)
  meta <- intersect(c("id", "gene_biotype", "Chr", "Start", "End", "Strand", "Length"), names(df))
  sample_cols <- setdiff(names(df), meta)
  if (!length(sample_cols)) stop("no sample columns in: ", path)
  clean <- clean_sample(sample_cols)
  names(df)[match(sample_cols, names(df))] <- clean
  lengths <- df$Length
  names(lengths) <- df$id
  counts <- df[, c("id", clean), drop = FALSE]
  list(counts = counts, lengths = lengths)
}

parse_count_matrix <- function(path, lengths) {
  df <- read_tsv(path)
  names(df)[1] <- "id"
  names(df)[-1] <- clean_sample(names(df)[-1])
  list(counts = df, lengths = lengths)
}

load_lengths <- function(path) {
  df <- read_tsv(path)
  if (ncol(df) < 2) stop("length table needs transcript and length columns: ", path)
  len <- as.numeric(df[[2]])
  names(len) <- as.character(df[[1]])
  len
}

pair_polysome <- function(samples, sheet, input_kw, poly_kw) {
  if (!is.null(sheet) && file.exists(sheet)) {
    meta <- read_tsv(sheet)
    names(meta) <- tolower(names(meta))
    need <- c("sample", "group", "fraction")
    if (!all(need %in% names(meta))) {
      stop("polysome sample sheet needs columns: sample, group, fraction")
    }
    meta$sample <- clean_sample(meta$sample)
    meta$fraction <- tolower(meta$fraction)
    inp <- meta[grepl(input_kw, meta$fraction, ignore.case = TRUE), , drop = FALSE]
    poly <- meta[grepl(poly_kw, meta$fraction, ignore.case = TRUE) | grepl("polysome", meta$fraction), , drop = FALSE]
    common <- intersect(inp$group, poly$group)
    if (!length(common)) stop("no Input/Polysome pairs in sample sheet groups")
    pairs <- do.call(rbind, lapply(common, function(g) {
      data.frame(
        pair = g,
        group = g,
        translation = poly$sample[poly$group == g][1],
        abundance = inp$sample[inp$group == g][1],
        stringsAsFactors = FALSE
      )
    }))
    return(pairs)
  }
  input_kw <- tolower(input_kw)
  poly_kw <- tolower(poly_kw)
  input_map <- list()
  poly_map <- list()
  for (s in samples) {
    sl <- tolower(s)
    if (grepl(input_kw, sl, fixed = TRUE)) {
      geno <- gsub(paste0("[_.-]?", input_kw, "[_.-]?"), "", sl)
      geno <- gsub("^[_-]+|[_-]+$", "", geno)
      if (!nzchar(geno)) geno <- s
      input_map[[geno]] <- s
    } else if (grepl(poly_kw, sl, fixed = TRUE)) {
      geno <- gsub(paste0("[_.-]?", poly_kw, "[_.-]?"), "", sl)
      geno <- gsub("^[_-]+|[_-]+$", "", geno)
      if (!nzchar(geno)) geno <- s
      poly_map[[geno]] <- s
    }
  }
  common <- intersect(names(input_map), names(poly_map))
  if (!length(common)) {
    stop("could not pair Input/Polysome samples with keywords '", input_kw, "' / '", poly_kw, "'")
  }
  do.call(rbind, lapply(sort(common), function(g) {
    data.frame(
      pair = g,
      group = g,
      translation = poly_map[[g]],
      abundance = input_map[[g]],
      stringsAsFactors = FALSE
    )
  }))
}

pair_ribo <- function(sheet) {
  if (is.null(sheet) || !file.exists(sheet)) {
    stop("Ribo-seq TE requires --sample-sheet with sample, group, assay, pair")
  }
  meta <- read_tsv(sheet)
  names(meta) <- tolower(names(meta))
  need <- c("sample", "group", "assay", "pair")
  if (!all(need %in% names(meta))) stop("Ribo sample sheet needs: sample, group, assay, pair")
  meta$sample <- clean_sample(meta$sample)
  meta$assay <- tolower(meta$assay)
  meta$assay[meta$assay == "qti"] <- "ribo"
  ribo <- meta[meta$assay == "ribo", , drop = FALSE]
  rna <- meta[meta$assay %in% c("rna", "input"), , drop = FALSE]
  common <- intersect(ribo$pair, rna$pair)
  if (!length(common)) stop("no ribo/rna pairs in sample sheet")
  do.call(rbind, lapply(common, function(p) {
    r1 <- ribo[ribo$pair == p, , drop = FALSE][1, ]
    r2 <- rna[rna$pair == p, , drop = FALSE][1, ]
    data.frame(
      pair = p,
      group = r1$group,
      translation = r1$sample,
      abundance = r2$sample,
      stringsAsFactors = FALSE
    )
  }))
}

pair_q <- function(pairs, te_ratio) {
  q <- rep(1, nrow(pairs))
  if (is.null(te_ratio)) {
    message("[INFO] no te_ratio table: manual TE uses classic TPM only (Q=1)")
    return(q)
  }
  for (i in seq_len(nrow(pairs))) {
    q_tr <- lookup_sf(te_ratio, pairs$translation[i])
    if (!is.finite(q_tr)) {
      stop(
        "te_ratio missing for translation sample ", pairs$translation[i],
        " (pair ", pairs$pair[i], ")"
      )
    }
    q[i] <- q_tr
    message(sprintf(
      "[INFO] pair %s: Q=%s on %s; Input/RNA remains 1",
      pairs$pair[i], format(q[i], digits = 6), pairs$translation[i]
    ))
  }
  q
}

mode <- get_arg("--mode", required = TRUE)
out_dir <- get_arg("--output-dir", required = TRUE)
spikein <- get_arg("--spikein-table", NULL)
sheet <- get_arg("--sample-sheet", NULL)
pseudo <- as.numeric(get_arg("--pseudocount", "1"))
min_tpm <- as.numeric(get_arg("--min-tpm", "1"))
te_ratio <- load_te_ratio(spikein)

if (mode == "polysome") {
  parsed <- parse_featurecounts(get_arg("--counts", required = TRUE))
  tpm <- tpm_from_counts(parsed$counts, parsed$lengths)
  sample_names <- names(tpm)[-1]
  pairs <- pair_polysome(
    sample_names, sheet,
    get_arg("--input-keyword", "input"),
    get_arg("--poly-keyword", "poly")
  )
  missing <- setdiff(c(pairs$translation, pairs$abundance), sample_names)
  if (length(missing)) stop("sample names not in count matrix: ", paste(missing, collapse = ", "))
  q <- pair_q(pairs, te_ratio)
  tpm_tr <- tpm[, c("id", pairs$translation), drop = FALSE]
  tpm_ab <- tpm[, c("id", pairs$abundance), drop = FALSE]
  names(tpm_tr)[-1] <- pairs$translation
  names(tpm_ab)[-1] <- pairs$abundance
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(tpm, file.path(out_dir, "TPM_matrix.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write_te(tpm$id, pairs, tpm_tr, tpm_ab, q, pseudo, min_tpm, out_dir, "gene")
} else if (mode == "ribo") {
  lengths <- load_lengths(get_arg("--lengths", required = TRUE))
  ribo <- parse_count_matrix(get_arg("--ribo-counts", required = TRUE), lengths)
  rna <- parse_count_matrix(get_arg("--rna-counts", required = TRUE), lengths)
  ribo_tpm <- tpm_from_counts(ribo$counts, lengths)
  rna_tpm <- tpm_from_counts(rna$counts, lengths)
  pairs <- pair_ribo(sheet)
  missing_r <- setdiff(pairs$translation, names(ribo_tpm)[-1])
  missing_n <- setdiff(pairs$abundance, names(rna_tpm)[-1])
  if (length(missing_r) || length(missing_n)) {
    stop(
      "sample sheet names missing from count matrices; ribo=",
      paste(missing_r, collapse = ","), " rna=", paste(missing_n, collapse = ",")
    )
  }
  q <- pair_q(pairs, te_ratio)
  common <- intersect(ribo_tpm$id, rna_tpm$id)
  ribo_tpm <- ribo_tpm[match(common, ribo_tpm$id), , drop = FALSE]
  rna_tpm <- rna_tpm[match(common, rna_tpm$id), , drop = FALSE]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ribo_tpm_q <- ribo_tpm
  for (i in seq_len(nrow(pairs))) {
    col <- pairs$translation[i]
    if (col %in% names(ribo_tpm_q)) ribo_tpm_q[[col]] <- ribo_tpm_q[[col]] * q[i]
  }
  utils::write.table(ribo_tpm_q, file.path(out_dir, "ribo_tpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(rna_tpm, file.path(out_dir, "rna_tpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  tpm_tr <- ribo_tpm[, c("id", pairs$translation), drop = FALSE]
  tpm_ab <- rna_tpm[, c("id", pairs$abundance), drop = FALSE]
  write_te(common, pairs, tpm_tr, tpm_ab, q, pseudo, min_tpm, out_dir, "transcript")
} else {
  stop("--mode must be polysome or ribo")
}

message("[INFO] manual TE written to ", out_dir)
