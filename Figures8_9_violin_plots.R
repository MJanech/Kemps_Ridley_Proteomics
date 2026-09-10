# ==============================================================================
# Figures8_9_violin_plots.R
#
# Purpose
#   Reproduce Figures 8 and 9 from the Kemp's ridley sea turtle plasma
#   proteomics manuscript.
#
# Input
#   Supplemental_Tables_1_through_8.xlsx
#
# Data source
#   - Figure values are original, untransformed weighted spectral counts (wSC)
#     from Table S3.
#   - BH-adjusted p-values are read directly from the paired Wilcoxon results:
#       Table S6 = OC2 admission vs release
#       Table S7 = OC4 admission vs release
#   - No abundance values are transformed or imputed for these figures.
#
# Figure structure
#   Figure 8: six proteins that increased during rehabilitation in at least one
#             oiling category.
#   Figure 9: eight proteins that decreased during rehabilitation in at least
#             one oiling category.
#
# Plot elements
#   - violin distributions for each of the four sampling groups
#   - individual raw wSC observations
#   - lines connecting paired admission/release samples from the same turtle
#   - dashed horizontal median lines
#   - BH-adjusted p-values from Tables S6 and S7
#
# Important
#   If a protein was excluded from a paired comparison by the manuscript's
#   comparison-specific zero-value filter, the plot is labeled "not tested"
#   for that comparison rather than assigning a p-value.
#
# Usage
#   Place this script and the supplemental workbook in the same directory:
#       source("Figures8_9_violin_plots.R")
#
#   Or provide the workbook path from the command line:
#       Rscript Figures8_9_violin_plots.R path/to/workbook.xlsx
# ==============================================================================

required_packages <- c("readxl", "ggplot2", "patchwork")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install these packages first: ",
    paste(missing_packages, collapse = ", "),
    "\nExample: install.packages(c(\"readxl\", \"ggplot2\", \"patchwork\"))"
  )
}

library(readxl)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# Input/output
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) {
  args[1]
} else {
  "Supplemental_Tables_1_through_8.xlsx"
}

if (!file.exists(input_file)) {
  stop(
    "Input workbook not found: ", input_file, "\n",
    "Place the workbook in the working directory or provide its path."
  )
}

output_dir <- "Figures8_9_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Helpers for reading the finalized Supporting Information workbook
# ------------------------------------------------------------------------------

normalize_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

find_column <- function(data, candidates) {
  normalized <- normalize_name(names(data))
  targets <- normalize_name(candidates)

  hit <- which(normalized %in% targets)

  if (length(hit) == 0) {
    return(NA_character_)
  }

  names(data)[hit[1]]
}

read_final_sheet <- function(sheet_name, skip_rows) {
  x <- read_excel(
    input_file,
    sheet = sheet_name,
    skip = skip_rows,
    .name_repair = "unique"
  )

  names(x) <- trimws(gsub("[\r\n]+", " ", names(x)))
  x
}

# Final workbook structure:
#   Table S3 header = Excel row 4 -> skip 3
#   Table S6 header = Excel row 4 -> skip 3
#   Table S7 header = Excel row 4 -> skip 3
s3 <- read_final_sheet("Table S3", 3)
s6 <- read_final_sheet("Table S6", 3)
s7 <- read_final_sheet("Table S7", 3)

s3_acc_col <- find_column(
  s3,
  c("Accession Number", "Accession.Number", "Accession")
)

s3_protein_col <- find_column(
  s3,
  c("Identified Proteins", "Identified.Proteins", "Protein")
)

s3_gene_col <- find_column(
  s3,
  c("Gene Name", "Gene.Name", "Gene")
)

s6_acc_col <- find_column(
  s6,
  c("Accession Number", "Accession.Number", "Accession")
)

s7_acc_col <- find_column(
  s7,
  c("Accession Number", "Accession.Number", "Accession")
)

s6_padj_col <- find_column(
  s6,
  c("padj_BH", "padj BH", "BH adjusted p-value")
)

s7_padj_col <- find_column(
  s7,
  c("padj_BH", "padj BH", "BH adjusted p-value")
)

required_detected <- c(
  s3_acc_col, s3_protein_col, s3_gene_col,
  s6_acc_col, s7_acc_col, s6_padj_col, s7_padj_col
)

if (any(is.na(required_detected))) {
  stop(
    "One or more required columns could not be detected in Tables S3, S6, or S7."
  )
}

# Trim accession strings because some finalized SI cells contain trailing spaces.
s3[[s3_acc_col]] <- trimws(as.character(s3[[s3_acc_col]]))
s6[[s6_acc_col]] <- trimws(as.character(s6[[s6_acc_col]]))
s7[[s7_acc_col]] <- trimws(as.character(s7[[s7_acc_col]]))

# ------------------------------------------------------------------------------
# Sample definitions
# ------------------------------------------------------------------------------

OC2_IDS <- c("114", "115", "116", "133", "137", "162", "192", "249", "278")

OC4_IDS <- c(
  "081", "087", "088", "098", "100", "107", "108", "110",
  "125", "193", "203", "204", "254", "265", "269"
)

GROUP_LEVELS <- c("OC2_0", "OC2_1", "OC4_0", "OC4_1")

sample_cols <- names(s3)[
  grepl("^Sample\\s*[0-9]{3}[AB]$", names(s3), ignore.case = TRUE)
]

if (length(sample_cols) != 48) {
  stop(
    "Expected 48 sample columns in Table S3 but found ",
    length(sample_cols), "."
  )
}

# ------------------------------------------------------------------------------
# Protein definitions and display precision
#
# The displayed number of decimal places matches the finalized manuscript
# figures. Exact BH-adjusted p-values are still written to the audit CSV.
# ------------------------------------------------------------------------------

fig8_spec <- data.frame(
  Panel = c("A", "B", "C", "D", "E", "F"),
  Accession = c(
    "XP_073178819.1",  # BGLAP
    "XP_073183428.1",  # CRP / pentraxin fusion protein-like
    "XP_073185772.1",  # COL1A2
    "XP_073189664.1",  # VMO1
    "XP_073166222.1",  # PLA2G4C
    "XP_073214170.1"   # ANGPTL3
  ),
  Plot_title = c(
    "osteocalcin",
    "pentraxin fusion protein-like",
    "collagen alpha-2(I) chain",
    "vitelline membrane outer layer protein 1 homolog",
    "Cluster of cytosolic phospholipase A2 gamma-like",
    "angiopoietin-related protein 3"
  ),
  OC2_digits = c(3, 2, 3, 2, 2, 2),
  OC4_digits = c(3, 2, 3, 3, 2, 3),
  stringsAsFactors = FALSE
)

fig9_spec <- data.frame(
  Panel = c("A", "B", "C", "D", "E", "F", "G", "H"),
  Accession = c(
    "XP_073170040.1",  # AK1
    "XP_073214567.1",  # FHL1
    "XP_073208099.1",  # PYGM
    "XP_073206612.1",  # MYL11
    "XP_073165987.1",  # AHCY
    "XP_073217254.1",  # PGK1
    "XP_073186293.1",  # MB
    "XP_073164022.1"   # GPI
  ),
  Plot_title = c(
    "Cluster of adenylate kinase isoenzyme 1 isoform X1",
    "four and a half LIM domains protein 1",
    "Cluster of glycogen phosphorylase, muscle form",
    "myosin regulatory light chain 11 isoform X1",
    "adenosylhomocysteinase",
    "phosphoglycerate kinase 1",
    "myoglobin isoform X1",
    "glucose-6-phosphate isomerase"
  ),
  OC2_digits = c(2, 2, 2, 2, 2, 2, 2, 2),
  OC4_digits = c(3, 3, 3, 2, 2, 3, 2, 2),
  stringsAsFactors = FALSE
)

all_specs <- rbind(
  transform(fig8_spec, Figure = "Figure 8"),
  transform(fig9_spec, Figure = "Figure 9")
)

# ------------------------------------------------------------------------------
# Retrieve exact BH-adjusted p-values from Tables S6 and S7
# ------------------------------------------------------------------------------

get_padj <- function(data, accession_col, padj_col, accession) {
  hit <- which(trimws(as.character(data[[accession_col]])) == accession)

  if (length(hit) == 0) {
    return(NA_real_)
  }

  if (length(hit) > 1) {
    warning(
      "Multiple rows found for accession ", accession,
      "; using the first row."
    )
  }

  suppressWarnings(as.numeric(data[[padj_col]][hit[1]]))
}

format_padj <- function(value, digits) {
  if (is.na(value)) {
    return(NA_character_)
  }

  formatC(
    value,
    format = "f",
    digits = digits
  )
}

stats_audit <- all_specs

stats_audit$OC2_padj_exact <- vapply(
  stats_audit$Accession,
  function(acc) get_padj(s6, s6_acc_col, s6_padj_col, acc),
  numeric(1)
)

stats_audit$OC4_padj_exact <- vapply(
  stats_audit$Accession,
  function(acc) get_padj(s7, s7_acc_col, s7_padj_col, acc),
  numeric(1)
)

stats_audit$OC2_label <- mapply(
  format_padj,
  stats_audit$OC2_padj_exact,
  stats_audit$OC2_digits,
  USE.NAMES = FALSE
)

stats_audit$OC4_label <- mapply(
  format_padj,
  stats_audit$OC4_padj_exact,
  stats_audit$OC4_digits,
  USE.NAMES = FALSE
)

write.csv(
  stats_audit,
  file.path(output_dir, "Figures8_9_statistics_audit.csv"),
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------------------------
# Figure colors
# ------------------------------------------------------------------------------

fill_colors <- c(
  "OC2_0" = "#F6C8D1",
  "OC2_1" = "#C9E2F2",
  "OC4_0" = "#F4E7B5",
  "OC4_1" = "#B8DED7"
)

point_colors <- c(
  "OC2_0" = "#D56F86",
  "OC2_1" = "#6FA9CF",
  "OC4_0" = "#C9A63E",
  "OC4_1" = "#399E95"
)

# ------------------------------------------------------------------------------
# Convert one Table S3 protein row to plotting data
# ------------------------------------------------------------------------------

protein_long <- function(accession) {
  hit <- which(s3[[s3_acc_col]] == accession)

  if (length(hit) == 0) {
    stop("Protein not found in Table S3: ", accession)
  }

  if (length(hit) > 1) {
    warning(
      "Multiple Table S3 rows found for accession ",
      accession, "; using the first row."
    )
  }

  hit <- hit[1]

  sample_name <- sub(
    "(?i)^Sample[^0-9]*",
    "",
    sample_cols,
    perl = TRUE
  )

  turtle_id <- sub("[AB]$", "", sample_name)

  occasion <- ifelse(
    grepl("A$", sample_name),
    "Admission",
    "Release"
  )

  oiling <- ifelse(
    turtle_id %in% OC2_IDS,
    "OC2",
    ifelse(
      turtle_id %in% OC4_IDS,
      "OC4",
      NA_character_
    )
  )

  group <- ifelse(
    oiling == "OC2" & occasion == "Admission", "OC2_0",
    ifelse(
      oiling == "OC2" & occasion == "Release", "OC2_1",
      ifelse(
        oiling == "OC4" & occasion == "Admission", "OC4_0",
        ifelse(
          oiling == "OC4" & occasion == "Release", "OC4_1",
          NA_character_
        )
      )
    )
  )

  wsc <- suppressWarnings(
    as.numeric(
      unlist(
        s3[hit, sample_cols],
        use.names = FALSE
      )
    )
  )

  if (anyNA(wsc)) {
    stop(
      "Unexpected missing/non-numeric wSC value for accession ",
      accession, "."
    )
  }

  data.frame(
    Accession = accession,
    Protein = as.character(s3[[s3_protein_col]][hit]),
    Gene = as.character(s3[[s3_gene_col]][hit]),
    Sample = sample_name,
    TurtleID = turtle_id,
    Occasion = occasion,
    OilingCategory = oiling,
    Group = factor(group, levels = GROUP_LEVELS),
    wSC = wsc,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Write the exact raw plotting data used for Figures 8 and 9
# ------------------------------------------------------------------------------

plotting_data <- do.call(
  rbind,
  lapply(seq_len(nrow(all_specs)), function(i) {
    d <- protein_long(all_specs$Accession[i])
    d$Figure <- all_specs$Figure[i]
    d$Panel <- all_specs$Panel[i]
    d$Plot_title <- all_specs$Plot_title[i]
    d
  })
)

write.csv(
  plotting_data,
  file.path(output_dir, "Figures8_9_raw_plotting_data.csv"),
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------------------------
# Single-panel plotting function
# ------------------------------------------------------------------------------

make_panel <- function(
  accession,
  title_text,
  oc2_label,
  oc4_label,
  show_x_title = FALSE,
  show_y_title = TRUE
) {

  d <- protein_long(accession)

  ymax <- max(d$wSC, na.rm = TRUE)

  if (!is.finite(ymax) || ymax <= 0) {
    ymax <- 1
  }

  top_limit <- ymax * 1.045
  annotation_y <- ymax * 0.915

  p <- ggplot(
    d,
    aes(x = Group, y = wSC)
  ) +
    geom_violin(
      aes(fill = Group),
      width = 0.70,
      scale = "area",
      trim = TRUE,
      colour = "grey25",
      linewidth = 0.40
    ) +
    geom_line(
      aes(group = interaction(OilingCategory, TurtleID)),
      colour = "grey45",
      linewidth = 0.36,
      alpha = 0.85
    ) +
    geom_point(
      aes(fill = Group),
      shape = 21,
      colour = "grey25",
      stroke = 0.28,
      size = 1.45,
      alpha = 0.95
    ) +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = 0.52,
      linewidth = 0.20,
      linetype = "dashed",
      colour = "grey25"
    ) +
    scale_fill_manual(
      values = fill_colors,
      guide = "none"
    ) +
    scale_x_discrete(
      drop = FALSE
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.00, 0.03))
    ) +
    coord_cartesian(
      ylim = c(0, top_limit),
      clip = "off"
    ) +
    labs(
      title = title_text,
      x = if (show_x_title) "Time of Collection" else NULL,
      y = if (show_y_title) "Weighted Spectral Count (wSC)" else NULL
    ) +
    theme_classic(
      base_family = "Arial",
      base_size = 9.5
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 9.4,
        hjust = 0.5,
        margin = margin(b = 1)
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 9.5,
        margin = margin(t = 1)
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 9.5,
        margin = margin(r = 1)
      ),
      axis.text.x = element_text(
        size = 8.3,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 8.3,
        colour = "black"
      ),
      axis.line = element_line(
        linewidth = 0.42,
        colour = "grey30"
      ),
      axis.ticks = element_line(
        linewidth = 0.38,
        colour = "grey30"
      ),
      plot.margin = margin(7, 5, 5, 7),
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      plot.background = element_rect(
        fill = "white",
        colour = NA
      )
    )

  # OC2 paired Wilcoxon BH FDR
  if (is.na(oc2_label)) {
    p <- p + annotate(
      "text",
      x = 1.42,
      y = annotation_y,
      label = "not tested",
      fontface = "italic",
      family = "Arial",
      size = 3.00
    )
  } else {
    p <- p + annotate(
      "text",
      x = 1.42,
      y = annotation_y,
      label = paste0("p[adj] == ", oc2_label),
      parse = TRUE,
      family = "Arial",
      size = 3.00
    )
  }

  # OC4 paired Wilcoxon BH FDR
  if (is.na(oc4_label)) {
    p <- p + annotate(
      "text",
      x = 3.42,
      y = annotation_y,
      label = "not tested",
      fontface = "italic",
      family = "Arial",
      size = 3.00
    )
  } else {
    p <- p + annotate(
      "text",
      x = 3.42,
      y = annotation_y,
      label = paste0("p[adj] == ", oc4_label),
      parse = TRUE,
      family = "Arial",
      size = 3.00
    )
  }

  p
}

# ------------------------------------------------------------------------------
# Build Figure 8
# ------------------------------------------------------------------------------

fig8_rows <- stats_audit[stats_audit$Figure == "Figure 8", ]

plots8 <- vector("list", nrow(fig8_rows))

for (i in seq_len(nrow(fig8_rows))) {
  plots8[[i]] <- make_panel(
    accession = fig8_rows$Accession[i],
    title_text = fig8_rows$Plot_title[i],
    oc2_label = fig8_rows$OC2_label[i],
    oc4_label = fig8_rows$OC4_label[i],
    show_x_title = i %in% c(5, 6),
    show_y_title = i %% 2 == 1
  )
}

figure8 <- wrap_plots(
  plots8,
  ncol = 2,
  byrow = TRUE
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 9.5
    ),
    plot.tag.position = c(0.005, 0.995)
  )

ggsave(
  file.path(output_dir, "Figure8_violin_plots.png"),
  figure8,
  width = 7.2,
  height = 9.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Figure8_violin_plots.pdf"),
  figure8,
  width = 7.2,
  height = 9.0,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

# ------------------------------------------------------------------------------
# Build Figure 9
# ------------------------------------------------------------------------------

fig9_rows <- stats_audit[stats_audit$Figure == "Figure 9", ]

plots9 <- vector("list", nrow(fig9_rows))

for (i in seq_len(nrow(fig9_rows))) {
  plots9[[i]] <- make_panel(
    accession = fig9_rows$Accession[i],
    title_text = fig9_rows$Plot_title[i],
    oc2_label = fig9_rows$OC2_label[i],
    oc4_label = fig9_rows$OC4_label[i],
    show_x_title = i %in% c(7, 8),
    show_y_title = i %% 2 == 1
  )
}

figure9 <- wrap_plots(
  plots9,
  ncol = 2,
  byrow = TRUE
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 9.5
    ),
    plot.tag.position = c(0.005, 0.995)
  )

ggsave(
  file.path(output_dir, "Figure9_violin_plots.png"),
  figure9,
  width = 7.2,
  height = 10.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Figure9_violin_plots.pdf"),
  figure9,
  width = 7.2,
  height = 10.8,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

# ------------------------------------------------------------------------------
# Reproducibility checks for the finalized figure statistics
# ------------------------------------------------------------------------------

# Figure 8
stopifnot(
  isTRUE(all.equal(
    stats_audit$OC2_padj_exact[
      stats_audit$Accession == "XP_073178819.1"
    ],
    0.0515919811320755,
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    stats_audit$OC4_padj_exact[
      stats_audit$Accession == "XP_073178819.1"
    ],
    0.000962690873579545,
    tolerance = 1e-12
  ))
)

# Figure 9 panel H (GPI) was retained in S6 but excluded from S7 by the
# comparison-specific zero-value filter. Therefore OC4 is correctly labeled
# "not tested".
stopifnot(
  !is.na(
    stats_audit$OC2_padj_exact[
      stats_audit$Accession == "XP_073164022.1"
    ]
  ),
  is.na(
    stats_audit$OC4_padj_exact[
      stats_audit$Accession == "XP_073164022.1"
    ]
  )
)

message("Figures 8 and 9 completed successfully.")
message("Data source: original, untransformed Table S3 wSC values.")
message("Statistics source: Tables S6 and S7 BH-adjusted paired Wilcoxon results.")
message("Outputs written to: ", normalizePath(output_dir))
