#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

parse_args <- function(values) {
  result <- list()
  for (i in seq(1L, length(values), by = 2L)) result[[sub("^--", "", values[[i]])]] <- values[[i + 1L]]
  result
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
lengths <- fread(args$lengths, sep = "\t")
offsets <- fread(args$offsets, sep = "\t")
metagene <- fread(args$metagene, sep = "\t")

save_both <- function(plot, stem, width = 7, height = 4.5) {
  ggsave(file.path(args$outdir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 180)
  ggsave(file.path(args$outdir, paste0(stem, ".pdf")), plot, width = width, height = height)
}

length_plot <- ggplot(lengths, aes(read_length, unique_alignments)) +
  geom_col(fill = "#2166AC") +
  labs(title = paste(args$library, "unique footprint lengths"), x = "Read length (nt)", y = "Unique alignments") +
  theme_bw(base_size = 11)
save_both(length_plot, "length_distribution")

offset_plot <- ggplot(offsets, aes(read_length, frame0_fraction, color = status)) +
  geom_hline(yintercept = 1 / 3, linetype = 2, color = "grey60") +
  geom_line(color = "grey50") + geom_point(size = 2) +
  geom_text(aes(label = recommended_offset), vjust = -0.7, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c(recommended = "#1B7837", insufficient_reads = "#B2182B")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = paste(args$library, "P-site offset recommendations"), x = "Read length (nt)",
       y = "Frame-0 fraction", color = "Status") + theme_bw(base_size = 11)
save_both(offset_plot, "offset_frame")

if (nrow(metagene)) {
  metagene_plot <- ggplot(metagene, aes(position, count, group = factor(read_length), color = factor(read_length))) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") + geom_line(linewidth = 0.5) +
    facet_wrap(~landmark, scales = "free_x") +
    labs(title = paste(args$library, "P-site metagene"), x = "Position relative to CDS landmark (nt)",
         y = "P-sites", color = "Length") + theme_bw(base_size = 11)
} else {
  metagene_plot <- ggplot() + annotate("text", x = 0, y = 0, label = "No recommended length/offset") +
    theme_void() + labs(title = paste(args$library, "P-site metagene"))
}
save_both(metagene_plot, "metagene", width = 8)

html <- c(
  "<!doctype html><html><head><meta charset='utf-8'><title>Ribo QC</title>",
  "<style>body{font-family:sans-serif;max-width:1100px;margin:2rem auto}img{max-width:95%;display:block;margin:1rem}</style></head><body>",
  sprintf("<h1>Ribo-seq QC: %s</h1>", args$library),
  "<p>Offsets are recommendations only. Freeze accepted lengths and offsets in the reviewed decision table.</p>",
  "<img src='length_distribution.png' alt='Length distribution'>",
  "<img src='offset_frame.png' alt='Offset and frame'>",
  "<img src='metagene.png' alt='P-site metagene'>",
  "</body></html>"
)
writeLines(html, file.path(args$outdir, "qc.html"))

