#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(patchwork)
  library(forcats)
  library(tibble)
  library(grid)
})

# =========================================================
# Input files
# =========================================================
len_file <- "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/tool_comparisions/HG00320-len-comp.tsv"
lev_file <- "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/tool_comparisions/HG00320-lev-comp.tsv"

n_loci_to_plot <- 3
selection_mode <- "top"   # "top" or "random"
set.seed(123)

out_pdf <- "supplement_examples_with_schematic_v5_1_large_labels.pdf"
out_png <- "supplement_examples_with_schematic_v5_1_large_labels.png"
out_csv <- "selected_informative_loci_with_schematic_v5_1_large_labels.csv"

# =========================================================
# Global readability settings
# =========================================================
BASE_SIZE <- 17
TITLE_SIZE <- 18
SUBTITLE_SIZE <- 14
AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 13
Y_AXIS_TEXT_SIZE <- 14
BAR_LABEL_SIZE <- 5.2
SCHEMATIC_TEXT_SIZE <- 5.2
SCHEMATIC_Y_TEXT_SIZE <- 14

# =========================================================
# Helpers
# =========================================================
parse_meta_lines <- function(file) {
  x <- readLines(file, warn = FALSE)
  
  meta_lines <- x[grepl("^##FILE_\\d+=<", x)]
  if (length(meta_lines) == 0) {
    stop("No ##FILE lines found in: ", file)
  }
  
  bind_rows(lapply(meta_lines, function(line) {
    tibble(
      file_idx = as.integer(str_match(line, "^##FILE_(\\d+)=")[, 2]),
      sample_name = str_trim(str_match(line, "Name:\\s*([^;>]+)")[, 2]),
      tool = str_trim(str_match(line, "Tool:\\s*([^;>]+)")[, 2]),
      start_offset = as.integer(str_match(line, "Start Offset:\\s*(-?\\d+)")[, 2]),
      end_offset = as.integer(str_match(line, "End Offset:\\s*(-?\\d+)")[, 2])
    )
  })) %>%
    arrange(file_idx)
}

read_comp_file <- function(file) {
  x <- readLines(file, warn = FALSE)
  header_idx <- which(!grepl("^##", x))[1]
  
  if (is.na(header_idx)) {
    stop("Could not find header line in: ", file)
  }
  
  read.delim(
    text = paste(x[header_idx:length(x)], collapse = "\n"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", ".", "")
  ) %>%
    mutate(
      CHROM = as.character(CHROM),
      START = as.numeric(START),
      END = as.numeric(END),
      MOT_LEN = as.numeric(MOT_LEN)
    )
}

make_key <- function(df) {
  df %>%
    mutate(locus_key = paste(CHROM, START, END, MOT_LEN, sep = "__"))
}

get_pair_columns <- function(df) {
  setdiff(names(df), c("CHROM", "START", "END", "MOT_LEN", "locus_key"))
}

count_nonzero_numeric <- function(df, cols) {
  apply(df[, cols, drop = FALSE], 1, function(x) {
    x <- suppressWarnings(as.numeric(x))
    sum(!is.na(x) & x != 0)
  })
}

count_nonmissing_numeric <- function(df, cols) {
  apply(df[, cols, drop = FALSE], 1, function(x) {
    x <- suppressWarnings(as.numeric(x))
    sum(!is.na(x))
  })
}

max_abs_numeric <- function(df, cols) {
  apply(df[, cols, drop = FALSE], 1, function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    max(abs(x))
  })
}

row_to_long_pairs <- function(row_df, pair_cols, value_name = "value") {
  tibble(
    pair = pair_cols,
    !!value_name := as.numeric(unlist(row_df[pair_cols]))
  ) %>%
    separate(pair, into = c("tool1", "tool2"), sep = "-", remove = FALSE, extra = "merge") %>%
    mutate(abs_value = abs(.data[[value_name]]))
}

shorten_tool_pair_labels <- function(x) {
  x %>%
    str_replace_all("Medaka Tandem", "Medaka") %>%
    str_replace_all("Straglr", "Straglr") %>%
    str_replace_all("LongTR", "LongTR") %>%
    str_replace_all("STRdust", "STRdust") %>%
    str_replace_all("STRkit", "STRkit") %>%
    str_replace_all("Atarva", "Atarva") %>%
    str_replace_all("vamos", "Vamos") %>%
    str_replace_all("ASM", "ASM")
}

# Left panels: first tool fixed as Straglr, second tool alphabetical
order_straglr_pairs <- function(df) {
  df %>%
    mutate(other_tool = ifelse(tool1 == "Straglr", tool2, tool1)) %>%
    arrange(other_tool) %>%
    mutate(pair = factor(pair, levels = rev(unique(pair))))
}

# Right panels: first tool alphabetical, then second tool alphabetical
order_lev_pairs <- function(df) {
  df %>%
    arrange(tool1, tool2) %>%
    mutate(pair = factor(pair, levels = rev(unique(pair))))
}

plot_pair_bar_panel <- function(long_df, title_text, xlab_text, use_absolute = FALSE) {
  df_plot <- long_df %>%
    filter(!is.na(value)) %>%
    mutate(
      plot_value = if (use_absolute) abs(value) else value,
      label = as.character(round(if (use_absolute) abs(value) else value, 0)),
      pair_label = shorten_tool_pair_labels(as.character(pair))
    )
  
  if (nrow(df_plot) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No data", size = 6) +
        labs(title = title_text, x = xlab_text, y = NULL) +
        theme_void(base_size = BASE_SIZE) +
        theme(plot.title = element_text(face = "bold", size = TITLE_SIZE))
    )
  }
  
  max_x <- max(df_plot$plot_value, na.rm = TRUE)
  x_pad <- max(3, max_x * 0.20)
  
  ggplot(df_plot, aes(x = plot_value, y = fct_reorder(pair_label, plot_value))) +
    geom_col(width = 0.84) +
    geom_text(
      aes(label = label),
      hjust = -0.06,
      size = BAR_LABEL_SIZE
    ) +
    labs(
      title = title_text,
      x = xlab_text,
      y = NULL
    ) +
    theme_bw(base_size = BASE_SIZE) +
    theme(
      plot.title = element_text(face = "bold", size = TITLE_SIZE, margin = margin(b = 12)),
      axis.title.x = element_text(size = AXIS_TITLE_SIZE, margin = margin(t = 10)),
      axis.text.x = element_text(size = AXIS_TEXT_SIZE),
      axis.text.y = element_text(size = Y_AXIS_TEXT_SIZE),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(t = 10, r = 50, b = 10, l = 14)
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.01, 0.26))
    ) +
    coord_cartesian(xlim = c(0, max_x + x_pad), clip = "off")
}

# =========================================================
# Schematic panels
# =========================================================
base_cols <- c("C" = "#F28E5B", "A" = "#69C3A5", "G" = "#8DA0CB")

make_repeat_df <- function(n_truth = 8, n_called = 9, motif = c("C", "A", "G"),
                           row_labels = c("Tool 1", "Tool 2")) {
  truth_seq <- rep(motif, length.out = n_truth * length(motif))
  called_seq <- rep(motif, length.out = n_called * length(motif))
  
  truth_df <- tibble(
    row = row_labels[1],
    idx = seq_along(truth_seq),
    base = truth_seq
  )
  
  called_df <- tibble(
    row = row_labels[2],
    idx = seq_along(called_seq),
    base = called_seq
  )
  
  bind_rows(truth_df, called_df) %>%
    mutate(
      x = idx,
      y = ifelse(row == row_labels[1], 2, 1)
    )
}

plot_schematic_length <- function() {
  df <- make_repeat_df(n_truth = 8, n_called = 9, row_labels = c("Tool 1", "Tool 2")) %>%
    mutate(is_excess = row == "Tool 2" & idx > 24)
  
  ggplot(df, aes(x = x, y = y)) +
    geom_tile(aes(fill = base), color = "white", width = 0.95, height = 0.72) +
    geom_text(aes(label = base), size = SCHEMATIC_TEXT_SIZE) +
    geom_tile(
      data = subset(df, is_excess),
      fill = NA, color = "red", linewidth = 1.3, width = 0.95, height = 0.72
    ) +
    scale_fill_manual(values = base_cols, guide = "none") +
    scale_y_continuous(breaks = c(1, 2), labels = c("Tool 2", "Tool 1")) +
    labs(
      title = "Absolute length difference",
      subtitle = "|called length - truth length| = |27 - 24| = 3 bp",
      x = NULL, y = NULL
    ) +
    coord_cartesian(xlim = c(0.5, 28.8), ylim = c(0.5, 2.5), clip = "off") +
    theme_void(base_size = BASE_SIZE) +
    theme(
      plot.title = element_text(face = "bold", size = TITLE_SIZE, margin = margin(b = 8)),
      plot.subtitle = element_text(size = SUBTITLE_SIZE, margin = margin(b = 10)),
      axis.text.y = element_text(size = SCHEMATIC_Y_TEXT_SIZE),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
}

plot_schematic_lev <- function() {
  truth <- c(rep(c("C", "A", "G"), 8))
  called <- c(rep(c("C", "A", "G"), 7), "C", "A", "A")
  
  truth_df <- tibble(row = "Tool 1", idx = seq_along(truth), base = truth)
  called_df <- tibble(row = "Tool 2", idx = seq_along(called), base = called)
  
  df <- bind_rows(truth_df, called_df) %>%
    mutate(
      y = ifelse(row == "Tool 1", 2, 1),
      is_error = row == "Tool 2" & idx == length(called)
    )
  
  ggplot(df, aes(x = idx, y = y)) +
    geom_tile(aes(fill = base), color = "white", width = 0.95, height = 0.72) +
    geom_text(aes(label = base), size = SCHEMATIC_TEXT_SIZE) +
    geom_tile(
      data = subset(df, is_error),
      fill = NA, color = "red", linewidth = 1.3, width = 0.95, height = 0.72
    ) +
    scale_fill_manual(values = base_cols, guide = "none") +
    scale_y_continuous(breaks = c(1, 2), labels = c("Tool 2", "Tool 1")) +
    labs(
      title = "Levenshtein distance",
      subtitle = "One substitution needed to convert Tool 2 to Tool 1: distance = 1",
      x = NULL, y = NULL
    ) +
    coord_cartesian(xlim = c(0.5, 24.8), ylim = c(0.5, 2.5), clip = "off") +
    theme_void(base_size = BASE_SIZE) +
    theme(
      plot.title = element_text(face = "bold", size = TITLE_SIZE, margin = margin(b = 8)),
      plot.subtitle = element_text(size = SUBTITLE_SIZE, margin = margin(b = 10)),
      axis.text.y = element_text(size = SCHEMATIC_Y_TEXT_SIZE),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
}

# =========================================================
# Read data
# =========================================================
len_meta <- parse_meta_lines(len_file)
lev_meta <- parse_meta_lines(lev_file)

len_df <- read_comp_file(len_file) %>% make_key() %>% distinct(locus_key, .keep_all = TRUE)
lev_df <- read_comp_file(lev_file) %>% make_key() %>% distinct(locus_key, .keep_all = TRUE)

len_pair_cols <- get_pair_columns(len_df)
lev_pair_cols <- get_pair_columns(lev_df)

lev_pair_cols_use <- lev_pair_cols
len_pair_cols_use <- len_pair_cols[
  str_detect(len_pair_cols, "^Straglr-") | str_detect(len_pair_cols, "-Straglr$")
]

# =========================================================
# Select informative loci
# =========================================================
common_keys <- intersect(len_df$locus_key, lev_df$locus_key)

len_info <- len_df %>%
  filter(locus_key %in% common_keys) %>%
  mutate(
    len_nonmissing = count_nonmissing_numeric(cur_data(), len_pair_cols_use),
    len_nonzero = count_nonzero_numeric(cur_data(), len_pair_cols_use),
    len_max_abs = max_abs_numeric(cur_data(), len_pair_cols_use),
    locus_size = END - START
  ) %>%
  select(locus_key, CHROM, START, END, MOT_LEN, locus_size, len_nonmissing, len_nonzero, len_max_abs)

lev_info <- lev_df %>%
  filter(locus_key %in% common_keys) %>%
  mutate(
    lev_nonmissing = count_nonmissing_numeric(cur_data(), lev_pair_cols_use),
    lev_nonzero = count_nonzero_numeric(cur_data(), lev_pair_cols_use),
    lev_max_abs = max_abs_numeric(cur_data(), lev_pair_cols_use)
  ) %>%
  select(locus_key, lev_nonmissing, lev_nonzero, lev_max_abs)

candidate_df <- len_info %>%
  inner_join(lev_info, by = "locus_key") %>%
  filter(
    len_nonmissing >= 1,
    lev_nonmissing >= 5,
    len_nonzero >= 1 | lev_nonzero >= 2,
    len_max_abs <= 500,
    lev_max_abs <= 500,
    len_max_abs >= 1 | lev_max_abs >= 2
  ) %>%
  mutate(
    score =
      (pmin(coalesce(len_max_abs, 0), 50) * 2) +
      (pmin(coalesce(lev_max_abs, 0), 50) * 2) +
      (len_nonzero + lev_nonzero) -
      coalesce(locus_size, 0) / 1000
  ) %>%
  arrange(desc(score), locus_size, CHROM, START) %>%
  distinct(locus_key, .keep_all = TRUE)

if (nrow(candidate_df) < n_loci_to_plot) {
  stop("Not enough informative non-extreme unique loci found after filtering.")
}

if (selection_mode == "random") {
  selected_keys <- sample(candidate_df$locus_key, n_loci_to_plot)
} else {
  selected_keys <- head(candidate_df$locus_key, n_loci_to_plot)
}

selected_len <- len_df %>%
  filter(locus_key %in% selected_keys) %>%
  distinct(locus_key, .keep_all = TRUE) %>%
  slice(match(selected_keys, locus_key))

selected_lev <- lev_df %>%
  filter(locus_key %in% selected_keys) %>%
  distinct(locus_key, .keep_all = TRUE) %>%
  slice(match(selected_keys, locus_key))

write.csv(
  candidate_df %>%
    filter(locus_key %in% selected_keys) %>%
    select(CHROM, START, END, MOT_LEN, locus_size, len_nonmissing, len_nonzero, len_max_abs,
           lev_nonmissing, lev_nonzero, lev_max_abs, score),
  out_csv,
  row.names = FALSE
)

# =========================================================
# Build example panels
# =========================================================
row_plots <- vector("list", length = n_loci_to_plot)

for (i in seq_len(n_loci_to_plot)) {
  len_row <- selected_len[i, , drop = FALSE]
  lev_row <- selected_lev[i, , drop = FALSE]
  
  len_long <- row_to_long_pairs(len_row, len_pair_cols_use, "value") %>%
    order_straglr_pairs()
  
  lev_long <- row_to_long_pairs(lev_row, lev_pair_cols_use, "value") %>%
    order_lev_pairs()
  
  locus_title <- paste0(
    len_row$CHROM, ":", format(len_row$START, big.mark = ","), "-", format(len_row$END, big.mark = ","),
    " | motif length = ", len_row$MOT_LEN, " bp"
  )
  
  p_len <- plot_pair_bar_panel(
    len_long,
    title_text = locus_title,
    xlab_text = "Absolute length difference (bp)",
    use_absolute = TRUE
  )
  
  p_lev <- plot_pair_bar_panel(
    lev_long,
    title_text = locus_title,
    xlab_text = "Levenshtein distance",
    use_absolute = FALSE
  )
  
  row_plots[[i]] <- p_len + p_lev + plot_layout(widths = c(1, 1.15))
}

schematic_row <- plot_schematic_length() + plot_schematic_lev() + plot_layout(widths = c(1, 1))
examples_block <- wrap_plots(row_plots, ncol = 1)

final_plot <- schematic_row / examples_block +
  plot_layout(heights = c(1.35, n_loci_to_plot * 1.8))

# =========================================================
# Save large high-resolution outputs
# =========================================================
pdf_height <- 5 + 5.8 * n_loci_to_plot
png_height <- 5 + 5.8 * n_loci_to_plot

ggsave(
  out_pdf,
  final_plot,
  width = 22.5,
  height = pdf_height,
  limitsize = FALSE,
  device = cairo_pdf
)

ggsave(
  out_png,
  final_plot,
  width = 22.5,
  height = png_height,
  dpi = 500,
  limitsize = FALSE
)

message("Done.")
message("Saved PDF: ", out_pdf)
message("Saved PNG: ", out_png)
message("Saved CSV: ", out_csv)

print(
  candidate_df %>%
    filter(locus_key %in% selected_keys) %>%
    select(CHROM, START, END, MOT_LEN, locus_size, len_max_abs, lev_max_abs, score)
)