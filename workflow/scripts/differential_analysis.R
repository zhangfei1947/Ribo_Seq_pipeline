#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
})

parse_args <- function(values) {
  result <- list()
  index <- 1L
  while (index <= length(values)) {
    key <- sub("^--", "", values[[index]])
    if (index == length(values)) stop("Missing value for --", key)
    result[[gsub("-", "_", key)]] <- values[[index + 1L]]
    index <- index + 2L
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("ribo", "rna_cds", "rna_exon", "libraries", "decisions", "contrasts",
              "study", "contrast", "outdir", "alpha", "lfc_threshold", "min_total_count")
missing_args <- setdiff(required, names(args))
if (length(missing_args)) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

alpha <- as.numeric(args$alpha)
lfc_threshold <- as.numeric(args$lfc_threshold)
min_total_count <- as.integer(args$min_total_count)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

read_matrix <- function(path) {
  table <- fread(path, sep = "\t", check.names = FALSE)
  stopifnot(all(c("gene_id", "gene_symbol") %in% names(table)))
  if (anyDuplicated(table$gene_id)) stop("Duplicate gene_id in ", path)
  table
}

ribo <- read_matrix(args$ribo)
rna_cds <- read_matrix(args$rna_cds)
rna_exon <- read_matrix(args$rna_exon)
if (!setequal(ribo$gene_id, rna_cds$gene_id)) stop("Ribo and RNA CDS matrices have different gene sets")
rna_cds <- rna_cds[match(ribo$gene_id, gene_id)]
libraries <- fread(args$libraries, sep = "\t")
decisions <- fread(args$decisions, sep = "\t")
contrast_table <- fread(args$contrasts, sep = "\t")
contrast_row <- contrast_table[study_id == args$study & contrast_id == args$contrast]
if (nrow(contrast_row) != 1L) stop("Expected exactly one requested contrast")
reference <- contrast_row$reference_condition[[1]]
test <- contrast_row$test_condition[[1]]

metadata <- merge(libraries, decisions, by = "library_id", all.x = TRUE)
metadata <- metadata[study_id == args$study & condition %in% c(reference, test)]
metadata[, condition := factor(condition, levels = c(reference, test))]
normalization_records <- list()

analysis_level <- function(meta) {
  counts <- table(meta$condition)
  minimum <- min(counts[c(reference, test)])
  if (minimum >= 3L) "formal" else if (minimum == 2L) "exploratory_low_replication" else "descriptive"
}

display_labels <- function(ids, symbols) {
  symbols[is.na(symbols) | symbols == ""] <- ids[is.na(symbols) | symbols == ""]
  duplicated_symbol <- duplicated(symbols) | duplicated(symbols, fromLast = TRUE)
  symbols[duplicated_symbol] <- paste0(symbols[duplicated_symbol], " [", ids[duplicated_symbol], "]")
  symbols
}

empty_result <- function(matrix, level, status) {
  data.table(
    gene_id = matrix$gene_id,
    gene_symbol = matrix$gene_symbol,
    baseMean = NA_real_, log2FoldChange = NA_real_, lfcSE = NA_real_, stat = NA_real_,
    pvalue = NA_real_, padj = NA_real_, analysis_level = level, significance = status
  )
}

descriptive_result <- function(matrix, meta, analysis_name) {
  sample_ids <- intersect(meta$library_id, setdiff(names(matrix), c("gene_id", "gene_symbol")))
  meta <- meta[match(sample_ids, library_id)]
  counts <- as.matrix(matrix[, ..sample_ids])
  storage.mode(counts) <- "integer"
  library_sizes <- colSums(counts)
  positive <- library_sizes > 0
  if (!all(positive)) stop("Zero-size count library: ", paste(sample_ids[!positive], collapse = ", "))
  factors <- library_sizes / exp(mean(log(library_sizes)))
  normalization_records[[length(normalization_records) + 1L]] <<- data.table(
    analysis = analysis_name, library_id = sample_ids, assay = meta$assay,
    normalization_method = "library_size_descriptive", size_factor = factors
  )
  normalized <- sweep(counts, 2L, factors, "/")
  ref_mean <- rowMeans(normalized[, meta$condition == reference, drop = FALSE])
  test_mean <- rowMeans(normalized[, meta$condition == test, drop = FALSE])
  data.table(
    gene_id = matrix$gene_id,
    gene_symbol = matrix$gene_symbol,
    baseMean = rowMeans(normalized),
    log2FoldChange = log2((test_mean + 0.5) / (ref_mean + 0.5)),
    lfcSE = NA_real_, stat = NA_real_, pvalue = NA_real_, padj = NA_real_,
    analysis_level = "descriptive", significance = "not_tested"
  )
}

valid_design <- function(formula, data) {
  matrix <- model.matrix(formula, data = data)
  qr(matrix)$rank == ncol(matrix)
}

run_deseq_robust <- function(dds) {
  tryCatch(
    DESeq(dds, quiet = TRUE),
    error = function(first_error) {
      message("Default dispersion fit failed; retrying fitType='mean': ", conditionMessage(first_error))
      tryCatch(
        DESeq(dds, fitType = "mean", quiet = TRUE),
        error = function(second_error) {
          message("Mean dispersion fit failed; using gene-wise dispersion estimates: ",
                  conditionMessage(second_error))
          fallback <- estimateSizeFactors(dds)
          fallback <- estimateDispersionsGeneEst(fallback, quiet = TRUE)
          dispersions(fallback) <- mcols(fallback)$dispGeneEst
          nbinomWaldTest(fallback, quiet = TRUE)
        }
      )
    }
  )
}

result_with_statistics <- function(raw, shrunk) {
  raw_table <- as.data.table(as.data.frame(raw), keep.rownames = "gene_id")
  effect_table <- as.data.table(as.data.frame(shrunk), keep.rownames = "gene_id")
  for (column in c("baseMean", "stat", "pvalue", "padj")) {
    if (!column %in% names(effect_table)) {
      effect_table[[column]] <- raw_table[[column]][match(effect_table$gene_id, raw_table$gene_id)]
    }
  }
  effect_table
}

run_condition_de <- function(matrix, assay_name, requested) {
  if (!requested) return(empty_result(matrix, "not_requested", "not_requested"))
  meta <- metadata[assay == assay_name & include_assay_analysis == "yes"]
  sample_ids <- intersect(meta$library_id, setdiff(names(matrix), c("gene_id", "gene_symbol")))
  meta <- meta[match(sample_ids, library_id)]
  if (!all(c(reference, test) %in% as.character(meta$condition))) {
    stop("Included ", assay_name, " libraries do not contain both contrast conditions")
  }
  level <- analysis_level(meta)
  if (level == "descriptive") return(descriptive_result(matrix, meta, paste0(assay_name, "_de")))
  counts <- as.matrix(matrix[, ..sample_ids])
  storage.mode(counts) <- "integer"
  rownames(counts) <- matrix$gene_id
  keep <- rowSums(counts) >= min_total_count
  coldata <- as.data.frame(meta)
  rownames(coldata) <- coldata$library_id
  design_formula <- ~ condition
  if ("batch" %in% names(coldata) && length(unique(coldata$batch)) > 1L) {
    candidate <- ~ batch + condition
    if (valid_design(candidate, coldata)) design_formula <- candidate
  }
  dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE], colData = coldata,
                                design = design_formula)
  dds <- run_deseq_robust(dds)
  normalization_records[[length(normalization_records) + 1L]] <<- data.table(
    analysis = paste0(assay_name, "_de"), library_id = names(sizeFactors(dds)), assay = assay_name,
    normalization_method = "DESeq2_median_ratio_within_assay", size_factor = as.numeric(sizeFactors(dds))
  )
  raw <- results(dds, contrast = c("condition", test, reference), alpha = alpha)
  shrunk <- tryCatch(
    lfcShrink(dds, contrast = c("condition", test, reference), res = raw, type = "normal", quiet = TRUE),
    error = function(e) raw
  )
  result <- result_with_statistics(raw, shrunk)
  result <- merge(matrix[, .(gene_id, gene_symbol)], result, by = "gene_id", all.x = TRUE, sort = FALSE)
  result[, `:=`(
    analysis_level = level,
    significance = fifelse(!is.na(padj) & padj < alpha & abs(log2FoldChange) >= lfc_threshold,
                           "significant", "not_significant")
  )]
  setcolorder(result, c("gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE",
                       "stat", "pvalue", "padj", "analysis_level", "significance"))
  result
}

descriptive_dte <- function(ribo_matrix, rna_matrix, paired_meta) {
  normalize <- function(counts, ids, assay_name) {
    sizes <- colSums(counts)
    factors <- sizes / exp(mean(log(sizes)))
    normalization_records[[length(normalization_records) + 1L]] <<- data.table(
      analysis = "dte", library_id = ids, assay = assay_name,
      normalization_method = "assay_separate_library_size_descriptive", size_factor = factors
    )
    sweep(counts, 2L, factors, "/")
  }
  ribo_ids <- paired_meta[assay == "ribo"]$library_id
  rna_ids <- paired_meta[assay == "rna"]$library_id
  ribo_counts <- normalize(as.matrix(ribo_matrix[, ..ribo_ids]), ribo_ids, "ribo")
  rna_counts <- normalize(as.matrix(rna_matrix[, ..rna_ids]), rna_ids, "rna")
  ribo_bio <- paired_meta[match(ribo_ids, library_id)]$bio_sample_id
  rna_bio <- paired_meta[match(rna_ids, library_id)]$bio_sample_id
  rna_counts <- rna_counts[, match(ribo_bio, rna_bio), drop = FALSE]
  te <- log2(ribo_counts + 0.5) - log2(rna_counts + 0.5)
  bio_condition <- paired_meta[match(ribo_bio, bio_sample_id)]$condition
  effect <- rowMeans(te[, bio_condition == test, drop = FALSE]) -
    rowMeans(te[, bio_condition == reference, drop = FALSE])
  data.table(
    gene_id = ribo_matrix$gene_id, gene_symbol = ribo_matrix$gene_symbol,
    baseMean = rowMeans(cbind(ribo_counts, rna_counts)), log2FoldChange = effect,
    lfcSE = NA_real_, stat = NA_real_, pvalue = NA_real_, padj = NA_real_,
    analysis_level = "descriptive", significance = "not_tested"
  )
}

run_dte <- function(requested) {
  if (!requested) return(empty_result(ribo, "not_requested", "not_requested"))
  eligible <- metadata[include_te == "yes"]
  pair_table <- dcast(eligible, bio_sample_id + condition ~ assay, value.var = "library_id")
  pair_table <- pair_table[!is.na(ribo) & !is.na(rna)]
  if (!all(c(reference, test) %in% as.character(pair_table$condition))) {
    stop("DTE requires at least one reviewed Ribo/RNA pair in both conditions")
  }
  paired_meta <- eligible[bio_sample_id %in% pair_table$bio_sample_id]
  level <- analysis_level(unique(paired_meta[, .(bio_sample_id, condition)]))
  if (level == "descriptive") return(descriptive_dte(ribo, rna_cds, paired_meta))

  ribo_ids <- pair_table$ribo
  rna_ids <- pair_table$rna
  ribo_counts <- as.matrix(ribo[, ..ribo_ids])
  rna_counts <- as.matrix(rna_cds[, ..rna_ids])
  colnames(ribo_counts) <- paste0(ribo_ids, "__ribo")
  colnames(rna_counts) <- paste0(rna_ids, "__rna")
  counts <- cbind(rna_counts, ribo_counts)
  storage.mode(counts) <- "integer"
  rownames(counts) <- ribo$gene_id
  coldata <- rbindlist(list(
    data.table(column = colnames(rna_counts), library_id = rna_ids, assay = "rna"),
    data.table(column = colnames(ribo_counts), library_id = ribo_ids, assay = "ribo")
  ))
  coldata <- merge(coldata, paired_meta[, .(library_id, bio_sample_id, condition, batch)], by = "library_id")
  coldata <- coldata[match(colnames(counts), column)]
  coldata[, bio_sample_id := factor(bio_sample_id)]
  coldata[, assay := factor(assay, levels = c("rna", "ribo"))]
  coldata[, condition := factor(condition, levels = c(reference, test))]
  # With one condition per biological sample, the condition main effect is absorbed
  # by sample blocking. A single indicator is the full-rank assay:condition contrast.
  coldata[, dte_interaction := factor(as.integer(assay == "ribo" & condition == test), levels = c(0, 1))]
  rownames(coldata) <- coldata$column
  keep <- rowSums(counts) >= min_total_count
  dds <- DESeqDataSetFromMatrix(
    countData = counts[keep, , drop = FALSE], colData = as.data.frame(coldata),
    design = ~ bio_sample_id + assay + dte_interaction
  )
  dds <- run_deseq_robust(dds)
  normalization_records[[length(normalization_records) + 1L]] <<- data.table(
    analysis = "dte", library_id = coldata$library_id, assay = as.character(coldata$assay),
    normalization_method = "DESeq2_median_ratio_combined_interaction", size_factor = as.numeric(sizeFactors(dds))
  )
  interaction_names <- grep("dte_interaction", resultsNames(dds), value = TRUE)
  if (length(interaction_names) != 1L) {
    stop("Could not identify unique assay-by-condition coefficient: ", paste(resultsNames(dds), collapse = ", "))
  }
  raw <- results(dds, name = interaction_names[[1]], alpha = alpha)
  shrunk <- tryCatch(
    lfcShrink(dds, coef = interaction_names[[1]], type = "apeglm", quiet = TRUE),
    error = function(e) tryCatch(
      lfcShrink(dds, coef = interaction_names[[1]], type = "normal", quiet = TRUE),
      error = function(e2) raw
    )
  )
  result <- result_with_statistics(raw, shrunk)
  result <- merge(ribo[, .(gene_id, gene_symbol)], result, by = "gene_id", all.x = TRUE, sort = FALSE)
  result[, `:=`(
    analysis_level = level,
    significance = fifelse(!is.na(padj) & padj < alpha & abs(log2FoldChange) >= lfc_threshold,
                           "significant", "not_significant")
  )]
  setcolorder(result, c("gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE",
                       "stat", "pvalue", "padj", "analysis_level", "significance"))
  result
}

ribo_result <- run_condition_de(ribo, "ribo", contrast_row$run_ribo_de[[1]] == "yes")
rna_result <- run_condition_de(rna_exon, "rna", contrast_row$run_rna_de[[1]] == "yes")
dte_result <- run_dte(contrast_row$run_dte[[1]] == "yes")

fwrite(ribo_result, file.path(args$outdir, "ribo_differential.tsv"), sep = "\t", na = "NA")
fwrite(rna_result, file.path(args$outdir, "rna_differential.tsv"), sep = "\t", na = "NA")
fwrite(dte_result, file.path(args$outdir, "dte_interaction.tsv"), sep = "\t", na = "NA")
if (length(normalization_records)) {
  fwrite(rbindlist(normalization_records, fill = TRUE), file.path(args$outdir, "normalization_factors.tsv"), sep = "\t")
}

integrated <- merge(
  rna_result[, .(gene_id, gene_symbol, rna_log2fc = log2FoldChange, rna_padj = padj,
                 rna_significance = significance)],
  ribo_result[, .(gene_id, ribo_log2fc = log2FoldChange, ribo_padj = padj,
                  ribo_significance = significance)],
  by = "gene_id", all = TRUE
)
integrated <- merge(
  integrated,
  dte_result[, .(gene_id, dte_log2fc = log2FoldChange, dte_padj = padj,
                 dte_significance = significance, analysis_level)],
  by = "gene_id", all = TRUE
)

integrated[, evidence_class := "no_evidence"]
inferential <- integrated$analysis_level %in% c("formal", "exploratory_low_replication")
rna_sig <- !is.na(integrated$rna_padj) & integrated$rna_padj < alpha & abs(integrated$rna_log2fc) >= lfc_threshold
ribo_sig <- !is.na(integrated$ribo_padj) & integrated$ribo_padj < alpha & abs(integrated$ribo_log2fc) >= lfc_threshold
dte_sig <- !is.na(integrated$dte_padj) & integrated$dte_padj < alpha & abs(integrated$dte_log2fc) >= lfc_threshold
same_direction <- sign(integrated$rna_log2fc) == sign(integrated$dte_log2fc)
expected_ribo_effect <- integrated$rna_log2fc + integrated$dte_log2fc

integrated[!inferential, evidence_class := "descriptive_no_inference"]
integrated[inferential & rna_sig & !dte_sig, evidence_class := "transcriptional"]
integrated[inferential & !rna_sig & dte_sig, evidence_class := "translation_only"]
integrated[inferential & rna_sig & dte_sig & same_direction, evidence_class := "intensified"]
integrated[inferential & rna_sig & dte_sig & !same_direction & ribo_sig, evidence_class := "buffered"]
integrated[inferential & rna_sig & dte_sig & !same_direction & !ribo_sig, evidence_class := "completely_buffered"]
integrated[inferential & !rna_sig & !dte_sig & ribo_sig, evidence_class := "rpf_only_not_dte"]
integrated[inferential & (is.na(rna_log2fc) | is.na(ribo_log2fc) | is.na(dte_log2fc)), evidence_class := "low_information"]
integrated[
  inferential & ribo_sig & (rna_sig | dte_sig) & !is.na(expected_ribo_effect) &
    abs(expected_ribo_effect) >= 0.25 & sign(ribo_log2fc) != sign(expected_ribo_effect),
  evidence_class := "complex"
]
fwrite(integrated, file.path(args$outdir, "integrated_results.tsv"), sep = "\t", na = "NA")

plot_volcano <- function(result, title, stem) {
  plot_data <- copy(result)
  plot_data[, display_label := display_labels(gene_id, gene_symbol)]
  plot_data[, minus_log10_padj := -log10(pmax(padj, .Machine$double.xmin))]
  plot_data[, significant := !is.na(padj) & padj < alpha & abs(log2FoldChange) >= lfc_threshold]
  labels <- plot_data[significant == TRUE][order(padj)][seq_len(min(10L, .N))]
  figure <- ggplot(plot_data, aes(log2FoldChange, minus_log10_padj, color = significant)) +
    geom_point(alpha = 0.55, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = 2, color = "grey55") +
    geom_hline(yintercept = -log10(alpha), linetype = 2, color = "grey55") +
    geom_text(data = labels, aes(label = display_label), check_overlap = TRUE, size = 2.6,
              vjust = -0.5, show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#7F7F7F"), na.value = "#D0D0D0") +
    labs(title = title, x = "log2 fold change", y = "-log10 adjusted p-value", color = "Significant") +
    theme_bw(base_size = 11)
  ggsave(file.path(args$outdir, paste0(stem, ".png")), figure, width = 7, height = 5, dpi = 180)
  ggsave(file.path(args$outdir, paste0(stem, ".pdf")), figure, width = 7, height = 5)
}

plot_volcano(ribo_result, paste(args$contrast, "Ribo-seq differential occupancy"), "ribo_volcano")
plot_volcano(rna_result, paste(args$contrast, "RNA-seq differential expression"), "rna_volcano")
plot_volcano(dte_result, paste(args$contrast, "differential translation efficiency"), "dte_volcano")

summary_table <- integrated[, .N, by = evidence_class][order(-N)]
html <- c(
  "<!doctype html><html><head><meta charset='utf-8'><title>Ribo-seq analysis</title>",
  "<style>body{font-family:sans-serif;max-width:1100px;margin:2rem auto}table{border-collapse:collapse}th,td{padding:.4rem .8rem;border:1px solid #ccc}img{max-width:32%}</style></head><body>",
  sprintf("<h1>%s — %s</h1>", args$study, args$contrast),
  sprintf("<p>Reference: %s; test: %s. Analysis level: %s.</p>", reference, test,
          paste(unique(na.omit(integrated$analysis_level)), collapse = ", ")),
  "<h2>Integrated evidence classes</h2><table><tr><th>Class</th><th>Genes</th></tr>",
  paste(sprintf("<tr><td>%s</td><td>%d</td></tr>", summary_table$evidence_class, summary_table$N), collapse = ""),
  "</table><h2>Effect plots</h2>",
  "<img src='rna_volcano.png' alt='RNA volcano'><img src='ribo_volcano.png' alt='Ribo volcano'><img src='dte_volcano.png' alt='DTE volcano'>",
  "<p>Normalization is relative. Without a validated external spike-in, these results do not establish an absolute global translation shift.</p>",
  "<p>Gene symbols are display labels; FBgn gene_id remains the analysis key. Unreplicated comparisons are descriptive and contain no p-values or significance calls.</p>",
  "</body></html>"
)
writeLines(html, file.path(args$outdir, "report.html"))
writeLines(c("complete", paste0("study=", args$study), paste0("contrast=", args$contrast)),
           file.path(args$outdir, "complete.txt"))
