# ============================================================
# FIXED VERSION — Reviewer corrections applied:
#   1. Figure 4A font size increased
#   2. Figure 4A x-axis label includes "(log10 scale)"
#   3. Bar count labels shifted right and centered with bars
#   4. Minor tick marks removed from panel A y-axis
#   5. Tool legend made taller / less crowded
#   6. Panel letters B/C shifted so they do not cover 100% labels
#   7. annotation_logticks removed to eliminate minor tick marks
# ============================================================

library(ggplot2)
library(patchwork)
library(scales)
library(ggpubr)
library(grid)
library(dplyr)

# -------------------------
# Tool palette + labels
# -------------------------
tableau_colors <- c(
  "#4E79A7", "#F28E2B", "#E15759",
  "#76B7B2", "#59A14F", "#EDC948", "#FF9DA7"
)

tools  <- c("longtr", "atarva", "medaka", "strkit", "vamos", "strdust", "straglr")
labels <- c("LongTR", "ATaRVa", "Medaka\nTandem", "STRkit", "vamos", "STRdust", "Straglr")

tool_levels <- tools
tool_labels <- setNames(labels, tools)
tool_colors <- setNames(tableau_colors, tools)

# -------------------------
# Themes
# -------------------------

# General theme for panels B and C
theme_like_second <- function(base_size = 14, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      axis.line  = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      
      legend.position   = "bottom",
      legend.title      = element_text(face = "bold"),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(0.9,  "cm"),
      
      plot.title  = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# Larger theme for panel A
theme_panel_A <- function(base_size = 16, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.title      = element_text(face = "bold", size = base_size),
      axis.text       = element_text(color = "black", size = base_size - 1),
      axis.text.y     = element_text(color = "black", size = base_size - 1),
      axis.line       = element_line(linewidth = 0.7, color = "black"),
      axis.ticks      = element_line(linewidth = 0.7, color = "black"),
      
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold", size = base_size - 1),
      legend.text        = element_text(size = base_size - 2),
      legend.key.height  = unit(0.95, "cm"),
      legend.key.width   = unit(1.0,  "cm"),
      legend.spacing.y   = unit(0.30, "cm"),
      legend.box.spacing = unit(0.30, "cm"),
      legend.margin      = margin(6, 6, 6, 6),
      
      plot.title  = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Force horizontal geom_text labels and adjust placement
force_text_angle0 <- function(p, label_size = 4.5) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomText")) {
      p$layers[[i]]$aes_params$angle <- 0
      p$layers[[i]]$aes_params$hjust <- -0.18
      p$layers[[i]]$aes_params$vjust <- 0.5
      p$layers[[i]]$aes_params$size  <- label_size
    }
  }
  p
}

# -------------------------
# Panel tags (A/B/C)
# -------------------------
tag_theme_A <- theme(
  plot.tag               = element_text(face = "bold", size = 18),
  plot.tag.position      = c(0.01, 0.99),
  plot.tag.justification = c(0, 1),
  plot.tag.margin        = margin(0, 4, 4, 0)
)

tag_theme_BC <- theme(
  plot.tag               = element_text(face = "bold", size = 18),
  plot.tag.position      = c(-0.06, 1.02),
  plot.tag.justification = c(0, 1),
  plot.tag.margin        = margin(0, 4, 4, 0)
)

# -------------------------
# Dynamic y-axis zoom helpers
# -------------------------
ymax_from_plot <- function(p) {
  b  <- ggplot_build(p)
  ys <- unlist(lapply(b$data, function(d) d$y))
  suppressWarnings(max(ys, na.rm = TRUE))
}

ymin_from_plot <- function(p) {
  b  <- ggplot_build(p)
  ys <- unlist(lapply(b$data, function(d) d$y))
  ys <- ys[is.finite(ys)]
  suppressWarnings(min(ys, na.rm = TRUE))
}

apply_y_percent_smart_dynamic <- function(p, pad = 0.01) {
  ymx <- ymax_from_plot(p)
  ymn <- ymin_from_plot(p)
  
  if (is.finite(ymx) && ymx > 1.5) {
    lo <- max(0, floor((ymn - pad * 100) / 5) * 5)
    hi <- 100
    p +
      scale_y_continuous(
        labels = function(x) paste0(round(x), "%"),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      coord_cartesian(ylim = c(lo, hi), clip = "off")
  } else {
    lo <- max(0, ymn - pad)
    hi <- 1
    p +
      scale_y_continuous(
        labels = percent_format(accuracy = 1),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      coord_cartesian(ylim = c(lo, hi), clip = "off")
  }
}

# ============================================================
# PANEL A — log10 bar chart
# ============================================================
strip_log_scales <- function(p) {
  if (length(p$scales$scales) == 0) return(p)
  keep <- list()
  for (s in p$scales$scales) {
    is_log <- FALSE
    if (!is.null(s$trans) && !is.null(s$trans$name)) {
      is_log <- grepl("log", s$trans$name, ignore.case = TRUE)
    }
    if (!is_log) keep <- c(keep, list(s))
  }
  p$scales$scales <- keep
  p
}

p_left_base <- p_good_counts |>
  force_text_angle0(label_size = 4.5) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "Count of loci (log10 scale)", fill = "Tool") +
  theme_panel_A(base_size = 16) +
  theme(
    axis.ticks.y.minor = element_blank(),
    axis.ticks.x.minor = element_blank(),
    plot.margin = margin(10, 45, 10, 10)
  )

p_left_base_safe <- strip_log_scales(p_left_base)

built_left <- ggplot_build(p_left_base_safe)
d0 <- built_left$data[[1]]

x_is_num <- is.numeric(d0$x) && any(is.finite(d0$x))
y_is_num <- is.numeric(d0$y) && any(is.finite(d0$y))

if (x_is_num && !y_is_num) {
  # x is the numeric (log) axis — annotation_logticks removed to suppress minor ticks
  p_left <- p_left_base_safe +
    scale_x_log10(
      breaks       = c(10, 100, 1000, 10000, 100000),
      labels       = scales::label_number(scale_cut = scales::cut_short_scale()),
      minor_breaks = NULL
    )
} else {
  # y is the numeric (log) axis — annotation_logticks removed to suppress minor ticks
  p_left <- p_left_base_safe +
    scale_y_log10(
      breaks       = c(10, 100, 1000, 10000, 100000),
      labels       = scales::label_number(scale_cut = scales::cut_short_scale()),
      minor_breaks = NULL
    )
}

p_left <- p_left +
  scale_fill_manual(
    values = tool_colors,
    breaks = tool_levels,
    labels = tool_labels[tool_levels],
    name   = "Tool",
    drop   = FALSE
  ) +
  labs(tag = "A") +
  tag_theme_A

# ============================================================
# PANEL B — Mendelian consistency by proband genotype
# ============================================================
p_good_alltools$data <- p_good_alltools$data %>%
  mutate(tool = factor(as.character(tool), levels = tool_levels))

p_rt <- p_good_alltools +
  labs(x = "Proband genotype", y = NULL) +
  scale_color_manual(
    values = tool_colors,
    breaks = tool_levels,
    labels = tool_labels[tool_levels],
    drop   = FALSE
  ) +
  guides(color = "none") +
  theme_like_second(base_size = 14)

p_rt <- apply_y_percent_smart_dynamic(p_rt) +
  labs(tag = "B") +
  tag_theme_BC +
  theme(plot.margin = margin(18, 14, 8, 18))

# ============================================================
# PANEL C — Mendelian consistency by motif length
# ============================================================
p_all$data <- p_all$data %>%
  mutate(tool = factor(as.character(tool), levels = tool_levels))

p_rb <- p_all +
  labs(x = "Motif length", y = NULL) +
  scale_color_manual(
    values = tool_colors,
    breaks = tool_levels,
    labels = tool_labels[tool_levels],
    drop   = FALSE
  ) +
  guides(color = "none") +
  theme_like_second(base_size = 14)

p_rb <- apply_y_percent_smart_dynamic(p_rb) +
  labs(tag = "C") +
  tag_theme_BC +
  theme(plot.margin = margin(18, 14, 8, 18))

# ============================================================
# Assemble: A | legend | (B / C)
# ============================================================
main_fig <- (p_left | guide_area() | (p_rt / p_rb)) +
  plot_layout(widths = c(1.15, 0.38, 1), guides = "collect") &
  theme(
    legend.position    = "right",
    legend.key.height  = unit(0.95, "cm"),
    legend.key.width   = unit(1.0, "cm"),
    legend.spacing.y   = unit(0.30, "cm"),
    legend.box.spacing = unit(0.30, "cm"),
    legend.margin      = margin(6, 6, 6, 6)
  )

ggsave(
  file.path(fig_dir, "MAIN_FIGURE_3panel_pubready_log10_midlegend.png"),
  main_fig, width = 13.5, height = 7.5, dpi = 300, bg = "white"
)

ggsave(
  file.path(fig_dir, "MAIN_FIGURE_3panel_pubready_log10_midlegend.pdf"),
  main_fig, width = 13.5, height = 7.5, bg = "white"
)

print(main_fig)