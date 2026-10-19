# ==============================================================================
# GSEA_Figures6_7.R
#
# Purpose
#   Generate manuscript Figures 6 and 7 from the collapsed Reactome GSEA results.
#
# Statistical context
#   - limma-trend is the PRIMARY protein-level differential-abundance analysis.
#   - Figures 6 and 7 display the separate nonparametric preranked GSEA results.
#   - No limma statistics are used in these pathway figures.
#
# Figure layout
#   - BH-FDR legend is outside the data panel.
#   - Figures 7A and 7B are sized for side-by-side use.
#   - Three BH-FDR legend circles are shown per figure.
#   - Common NES x-axis limits are used across all pathway figures.
#   - PNG files are written at 600 dpi.
# ==============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2 first with: install.packages('ggplot2')")
}
if (!requireNamespace("grid", quietly = TRUE)) {
  stop("The grid package is required.")
}

library(ggplot2)
library(grid)

input_dir <- "NP_ranked_fgsea_collapsed_main_pathways"
out_dir <- "Revised_Figures6_7_NP_GSEA"

if (!dir.exists(input_dir)) {
  stop(paste0("Cannot find folder: ", input_dir,
              "\nPlace this script in the same working directory as that folder."))
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

find_input_file <- function(filename) {
  direct <- file.path(input_dir, filename)
  if (file.exists(direct)) return(direct)
  hits <- list.files(input_dir, pattern = paste0("^", filename, "$"), recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop(paste("Could not find", filename, "inside", input_dir))
  hits[1]
}

wrap_label <- function(x, width = 34) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

clean_pathway_label <- function(x) {
  raw <- gsub("^REACTOME_", "", x, ignore.case = TRUE)
  raw <- gsub("_", " ", raw)
  key <- tolower(raw)
  key <- gsub("[/\\-]", " ", key)
  key <- gsub("[^a-z0-9 ]", " ", key)
  key <- gsub("\\s+", " ", trimws(key))

  y <- tools::toTitleCase(tolower(raw))
  y[grepl("formation of fibrin clot", key)] <- "Fibrin Clot / Clotting Cascade"
  y[grepl("metabolism of carbohydrates and carbohydrate derivatives", key)] <- "Carbohydrate Metabolism"
  y[grepl("platelet aggregation", key) & grepl("plug formation", key)] <- "Platelet Aggregation / Plug Formation"
  y[grepl("creation of c4 and c2 activators", key)] <- "Creation of C4/C2 Activators"
  y[grepl("post translational protein modification", key)] <- "Post-Translational Protein Modification"
  y[grepl("l1cam interactions", key)] <- "L1CAM Interactions"
  y[grepl("signaling by nuclear receptors", key)] <- "Signaling by Nuclear Receptors"
  y[grepl("adaptive immune system", key)] <- "Adaptive Immune System"
  y[grepl("syndecan interactions", key)] <- "Syndecan Interactions"
  y[grepl("complement cascade", key)] <- "Complement Cascade"
  y[grepl("glucose metabolism", key)] <- "Glucose Metabolism"
  wrap_label(y, width = 34)
}

read_main <- function(comp) {
  f <- find_input_file(paste0(comp, "_collapsed_main_pathways.csv"))
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(d) == 0) return(d)
  req <- c("pathway", "padj", "NES", "leadingEdge")
  miss <- setdiff(req, names(d))
  if (length(miss) > 0) stop(paste("Missing required columns in", basename(f), ":", paste(miss, collapse = ", ")))
  d$Plot_label <- if ("Pathway_label" %in% names(d)) clean_pathway_label(d$Pathway_label) else clean_pathway_label(d$pathway)
  d$minus_log10_FDR <- -log10(pmax(d$padj, .Machine$double.xmin))
  d$Comparison <- comp
  d
}

add_direction <- function(d, comp) {
  if (nrow(d) == 0) return(d)
  if (comp == "S4") d$Biological_direction <- ifelse(d$NES > 0, "enriched in OC4", "enriched in OC2")
  if (comp %in% c("S6", "S7")) d$Biological_direction <- ifelse(d$NES > 0, "enriched at release", "enriched at admission")
  d
}

compute_size_info <- function(x, range_sizes = c(4.5, 9)) {
  x <- sort(as.numeric(x[is.finite(x)]))
  xmin <- min(x); xmax <- max(x)
  if (length(x) == 0) {
    breaks <- c(1, 2, 3)
  } else if (abs(xmax - xmin) < 1e-12) {
    delta <- max(0.05, 0.1 * max(xmax, 1))
    breaks <- c(max(0, xmin - delta), xmin, xmin + delta)
  } else {
    breaks <- c(xmin, (xmin + xmax) / 2, xmax)
  }
  labels <- format(signif(breaks, 3), trim = TRUE, scientific = FALSE)
  if (length(unique(labels)) < 3) labels <- format(signif(breaks, 4), trim = TRUE, scientific = FALSE)
  if (length(unique(labels)) < 3) labels <- formatC(breaks, digits = 3, format = "f")
  legend_sizes <- c(range_sizes[1], mean(range_sizes), range_sizes[2])
  list(breaks = breaks, labels = labels, legend_sizes = legend_sizes)
}

s4 <- add_direction(read_main("S4"), "S4")
s6 <- add_direction(read_main("S6"), "S6")
s7 <- add_direction(read_main("S7"), "S7")

all_nes <- c(s4$NES, s6$NES, s7$NES)
common_limit <- max(abs(all_nes), na.rm = TRUE)
common_limit <- ceiling((common_limit + 0.15) * 10) / 10
common_xlim <- c(-common_limit, common_limit)

if (.Platform$OS.type == "windows") {
  try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
}

base_family <- "Arial"
base_size <- 13
point_range <- c(4.5, 9)
figure6_width <- 10.2
figure7_width <- 7.3

make_plot <- function(d, title_text, subtitle_text, outfile_stub, fig_width = 10, y_wrap = 34) {
  if (nrow(d) == 0) return(invisible(NULL))

  d <- d[order(d$NES), , drop = FALSE]
  d$Plot_label <- factor(d$Plot_label, levels = d$Plot_label)
  size_info <- compute_size_info(d$minus_log10_FDR, range_sizes = point_range)

  p <- ggplot(d, aes(x = NES, y = Plot_label)) +
    geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed") +
    geom_segment(aes(x = 0, xend = NES, yend = Plot_label), linewidth = 0.8, colour = "grey40") +
    geom_point(aes(size = minus_log10_FDR, colour = Biological_direction)) +
    scale_x_continuous(limits = common_xlim) +
    scale_size_continuous(
      name = expression(-log[10]~"BH FDR"),
      range = point_range,
      limits = range(d$minus_log10_FDR, na.rm = TRUE),
      breaks = size_info$breaks,
      labels = size_info$labels,
      guide = guide_legend(
        order = 1,
        title.position = "top",
        override.aes = list(size = size_info$legend_sizes)
      )
    ) +
    guides(colour = "none") +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.position = "right",
      legend.justification = c(0, 0.5),
      legend.box.spacing = unit(0.10, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.key = element_blank(),
      legend.background = element_blank(),
      plot.margin = margin(12, 4, 10, 12)
    )

  fig_height <- max(4.8, 0.62 * nrow(d) + 2.5)

  ggsave(file.path(out_dir, paste0(outfile_stub, ".png")), p,
         width = fig_width, height = fig_height,
         units = "in", dpi = 600, bg = "white")
  ggsave(file.path(out_dir, paste0(outfile_stub, ".pdf")), p,
         width = fig_width, height = fig_height,
         units = "in", bg = "white", device = grDevices::cairo_pdf)

  cat(outfile_stub, " width = ", fig_width, " in; legend labels: ", paste(size_info$labels, collapse = ", "), "\n", sep = "")
  invisible(p)
}

make_plot(
  s4,
  "OC2 vs OC4 at Admission",
  "Positive NES = enriched in OC4; negative NES = enriched in OC2",
  "Figure6_S4_NP_GSEA_main_pathways",
  fig_width = figure6_width
)

make_plot(
  s6,
  "OC2 Admission vs Release",
  "Positive NES = enriched at release; negative NES = enriched at admission",
  "Figure7A_S6_OC2_NP_GSEA_main_pathways",
  fig_width = figure7_width
)

make_plot(
  s7,
  "OC4 Admission vs Release",
  "Positive NES = enriched at release; negative NES = enriched at admission",
  "Figure7B_S7_OC4_NP_GSEA_main_pathways",
  fig_width = figure7_width
)

plot_cols <- c("Comparison", "pathway", "Plot_label", "NES", "pval", "padj",
               "minus_log10_FDR", "Biological_direction", "leadingEdge")
plot_data <- rbind(
  s4[, intersect(plot_cols, names(s4)), drop = FALSE],
  s6[, intersect(plot_cols, names(s6)), drop = FALSE],
  s7[, intersect(plot_cols, names(s7)), drop = FALSE]
)
write.csv(plot_data, file.path(out_dir, "Revised_Figures6_7_exact_plotting_data.csv"), row.names = FALSE)

cat("\nUPDATED REVISED FIGURE 6/7 PLOTTING COMPLETE\n")
cat("==============================================\n\n")
cat("Common x-axis limits applied to all figures: ", common_xlim[1], " to ", common_xlim[2], "\n", sep = "")
cat("Figure 6 kept wide; Figures 7A and 7B narrowed for side-by-side use.\n")
cat("BH-FDR legend placed outside the panel and close to the plot body.\n")
cat("Larger x-axis and legend fonts retained.\n")
cat("Output folder: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n\n", sep = "")
cat("PNG files are 600 dpi; matching PDF files were also written.\n")
cat("No limma statistics were used in these pathway figures.\n")

