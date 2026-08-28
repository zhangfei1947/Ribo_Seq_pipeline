suppressPackageStartupMessages(library(data.table))

root <- normalizePath(file.path(getwd()))
tmp <- tempfile("riboseq-smoke-")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

genes <- sprintf("FBgn%07d", 1:30)
symbols <- paste0("Gene", 1:30)
bios <- c("C1", "C2", "T1", "T2")
conditions <- c("control", "control", "treated", "treated")
ribo_ids <- paste0("RIBO_", bios)
rna_ids <- paste0("RNA_", bios)

set.seed(11)
base <- matrix(rnbinom(30 * 4, mu = 120, size = 20), nrow = 30)
rna_counts <- base
ribo_counts <- base
rna_counts[1:5, 3:4] <- rna_counts[1:5, 3:4] * 2L
ribo_counts[6:10, 3:4] <- ribo_counts[6:10, 3:4] * 3L

write_matrix <- function(path, ids, counts) {
  table <- data.table(gene_id = genes, gene_symbol = symbols)
  for (i in seq_along(ids)) table[[ids[[i]]]] <- as.integer(counts[, i])
  fwrite(table, path, sep = "\t")
}
write_matrix(file.path(tmp, "ribo.tsv"), ribo_ids, ribo_counts)
write_matrix(file.path(tmp, "rna_cds.tsv"), rna_ids, rna_counts)
write_matrix(file.path(tmp, "rna_exon.tsv"), rna_ids, rna_counts)

libraries <- rbindlist(lapply(seq_along(bios), function(i) {
  rbind(
    data.table(library_id = ribo_ids[[i]], study_id = "study", bio_sample_id = bios[[i]], assay = "ribo", condition = conditions[[i]], layout = "SE", batch = "b1", strandedness = "forward", umi_pattern = "none", ribo_protocol_group = "p1", adapter_3p = "none", run_ids = ribo_ids[[i]], n_runs = 1),
    data.table(library_id = rna_ids[[i]], study_id = "study", bio_sample_id = bios[[i]], assay = "rna", condition = conditions[[i]], layout = "SE", batch = "b1", strandedness = "forward", umi_pattern = "none", ribo_protocol_group = "", adapter_3p = "none", run_ids = rna_ids[[i]], n_runs = 1)
  )
}))
fwrite(libraries, file.path(tmp, "libraries.tsv"), sep = "\t")
decisions <- libraries[, .(library_id, include_assay_analysis = "yes", include_te = "yes")]
fwrite(decisions, file.path(tmp, "decisions.tsv"), sep = "\t")
contrasts <- data.table(study_id = "study", contrast_id = "treated_vs_control", test_condition = "treated", reference_condition = "control", run_ribo_de = "yes", run_rna_de = "yes", run_dte = "yes")
fwrite(contrasts, file.path(tmp, "contrasts.tsv"), sep = "\t")

outdir <- file.path(tmp, "out")
status <- system2(
  "Rscript",
  shQuote(c(
    file.path(root, "workflow/scripts/differential_analysis.R"),
    "--ribo", file.path(tmp, "ribo.tsv"), "--rna-cds", file.path(tmp, "rna_cds.tsv"),
    "--rna-exon", file.path(tmp, "rna_exon.tsv"), "--libraries", file.path(tmp, "libraries.tsv"),
    "--decisions", file.path(tmp, "decisions.tsv"), "--contrasts", file.path(tmp, "contrasts.tsv"),
    "--study", "study", "--contrast", "treated_vs_control", "--outdir", outdir,
    "--alpha", "0.05", "--lfc-threshold", "1", "--min-total-count", "10"
  ))
)
stopifnot(status == 0L)
stopifnot(file.exists(file.path(outdir, "complete.txt")))
stopifnot(file.exists(file.path(outdir, "normalization_factors.tsv")))
stopifnot(all(c("gene_id", "gene_symbol", "dte_log2fc") %in% names(fread(file.path(outdir, "integrated_results.tsv")))))
cat("Differential smoke test passed\n")
