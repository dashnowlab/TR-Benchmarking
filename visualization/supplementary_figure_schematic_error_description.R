library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# ----------------------------
# User-adjustable parameters
# ----------------------------
motif <- "CAG"
motif_len <- nchar(motif)
true_repeats <- 8

example_df <- tibble(
  category = c("Match", "Off by 1 bp", "Off by motif", "Mismatch"),
  call_bp_diff = c(
    0,
    1,
    motif_len,
    motif_len + 2
  )
)

# ----------------------------
# Helpers
# ----------------------------
make_true_seq <- function(motif, n_rep) {
  paste0(rep(motif, n_rep), collapse = "")
}

make_called_seq <- function(motif, n_rep, bp_diff) {
  base_seq <- make_true_seq(motif, n_rep)
  
  if (bp_diff == 0) {
    return(base_seq)
  } else if (bp_diff > 0) {
    extra <- substr(paste0(rep(motif, 10), collapse = ""), 1, bp_diff)
    return(paste0(base_seq, extra))
  } else {
    keep_n <- max(0, nchar(base_seq) + bp_diff)
    return(substr(base_seq, 1, keep_n))
  }
}

seq_to_tiles <- function(seq_string, row_name, panel_name, y_value) {
  chars <- strsplit(seq_string, "")[[1]]
  tibble(
    x = seq_along(chars),
    base = chars,
    row = row_name,
    panel = panel_name,
    y = y_value
  )
}

# ----------------------------
# Build example sequences
# ----------------------------
true_seq <- make_true_seq(motif, true_repeats)

plot_df <- example_df %>%
  mutate(
    category = factor(
      category,
      levels = c("Match", "Off by 1 bp", "Off by motif", "Mismatch")
    ),
    truth_seq  = true_seq,
    called_seq = map_chr(call_bp_diff, ~ make_called_seq(motif, true_repeats, .x)),
    truth_len  = nchar(truth_seq),
    called_len = nchar(called_seq),
    abs_diff   = abs(called_len - truth_len),
    label = case_when(
      category == "Match" ~ "No error",
      category == "Off by 1 bp" ~ "Error = 1 bp",
      category == "Off by motif" ~ paste0("Error = ", abs_diff, " bp (<= motif length)"),
      category == "Mismatch" ~ paste0("Mismatch = ", abs_diff, " bp (> motif length)")
    )
  )

# ----------------------------
# Sequence tiles
# ----------------------------
tile_df <- bind_rows(
  map2_dfr(plot_df$truth_seq, plot_df$category,
           ~ seq_to_tiles(.x, "Truth", .y, 2)),
  map2_dfr(plot_df$called_seq, plot_df$category,
           ~ seq_to_tiles(.x, "Called", .y, 1))
)

# ----------------------------
# Use ONE global panel width
# so Truth has same visual size in all panels
# ----------------------------
global_max_len <- max(plot_df$called_len, plot_df$truth_len)

panel_info <- plot_df %>%
  transmute(
    panel = category,
    truth_len,
    called_len,
    panel_len = global_max_len,
    label
  )

pad_df <- panel_info %>%
  rowwise() %>%
  do(
    tibble(
      panel = .$panel,
      x = seq_len(.$panel_len)
    )
  ) %>%
  ungroup()

truth_pad <- pad_df %>%
  left_join(panel_info, by = "panel") %>%
  transmute(
    panel, x,
    row = "Truth",
    y = 2,
    exists = x <= truth_len
  )

called_pad <- pad_df %>%
  left_join(panel_info, by = "panel") %>%
  transmute(
    panel, x,
    row = "Called",
    y = 1,
    exists = x <= called_len
  )

pad_long <- bind_rows(truth_pad, called_pad) %>%
  left_join(tile_df, by = c("panel", "x", "row", "y")) %>%
  left_join(panel_info, by = "panel") %>%
  mutate(
    fill_group = case_when(
      !exists ~ "Empty",
      base == "A" ~ "A",
      base == "C" ~ "C",
      base == "G" ~ "G",
      base == "T" ~ "T",
      TRUE ~ "Other"
    ),
    display_base = ifelse(exists, base, ""),
    overhang = case_when(
      row == "Called" & exists & x > truth_len & called_len > truth_len ~ TRUE,
      row == "Truth"  & exists & x > called_len & truth_len > called_len ~ TRUE,
      TRUE ~ FALSE
    )
  )

ann_df <- panel_info %>%
  transmute(
    panel,
    x = (panel_len + 1) / 2,
    y = 2.9,
    label
  )

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(pad_long, aes(x = x, y = y)) +
  geom_tile(
    data = subset(pad_long, exists),
    aes(fill = fill_group),
    color = NA,
    width = 0.95,
    height = 0.75
  ) +
  geom_tile(
    data = subset(pad_long, overhang),
    fill = NA,
    color = "red3",
    linewidth = 1.2,
    width = 0.95,
    height = 0.75
  ) +
  geom_text(
    data = subset(pad_long, exists),
    aes(label = display_base),
    size = 3.2,
    family = "mono"
  ) +
  geom_text(
    data = ann_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  facet_wrap(~panel, ncol = 1, drop = FALSE) +
  scale_y_continuous(
    breaks = c(1, 2),
    labels = c("Called", "Truth"),
    expand = expansion(mult = c(0.05, 0.18))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_fill_manual(
    values = c(
      A = "#66c2a5",
      C = "#fc8d62",
      G = "#8da0cb",
      T = "#e78ac3",
      Empty = "white",
      Other = "grey80"
    )
  ) +
  labs(
    x = "Allele sequence",
    y = NULL,
    title = "Schematic definition of error categories",
    subtitle = paste0(
      "Example locus with motif = ", motif,
      " (motif length = ", motif_len, " bp). Red outline marks the excess."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p)

ggsave("supplementary_mendelian_error_schematic_fixed_truth.png", p, width = 8, height = 9, dpi = 300)
ggsave("supplementary_mendelian_error_schematic_fixed_truth.pdf", p, width = 8, height = 9)