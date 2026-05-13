# ============================================================
# FULL R SCRIPT (COPY/PASTE)
# Resource usage boxplots (Memory + Real-time) with:
#   - your colors + labels
#   - tools ordered alphabetically (by tool KEY)
#   - Straglr split into its own mini-panel (different scale)
#   - panel tags A and B (A=Memory, B=Real-time)
#   - ONE shared legend at bottom (horizontal)
#   - STACKED: A above B
#
# Expected columns (case-insensitive):
#   process, peak_rss, realtime
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# -------------------------
# USER INPUT
# -------------------------
infile <- "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/Supplementary_Table3_resources.tsv"
sep_is_tab <- TRUE
mem_col <- "peak_rss"   # <<< CHANGED (was rss)

# -------------------------
# Your palette + labels
# -------------------------
tableau_colors  <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                     "#59A14F", "#EDC948", "#FF9DA7")
tools_given  <- c("longtr", "atarva", "medaka", "strkit", "vamos", "strdust", "straglr")
labels_given <- c("LongTR", "ATaRVa", "Medaka\nTandem", "STRkit", "vamos", "STRdust", "Straglr")

tools_alpha <- sort(tools_given)
tools_left  <- sort(setdiff(tools_alpha, "straglr"))

label_map <- setNames(labels_given, tools_given)[tools_alpha]
color_map <- setNames(tableau_colors, tools_given)[tools_alpha]

# -------------------------
# Helpers
# -------------------------
parse_mem_to_gb <- function(x) {
  x <- str_trim(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  val  <- suppressWarnings(as.numeric(str_extract(x, "[0-9]+\\.?[0-9]*")))
  unit <- str_to_upper(str_extract(x, "(KB|MB|GB|TB)"))
  dplyr::case_when(
    is.na(val)   ~ NA_real_,
    unit == "KB" ~ val / (1024^2),
    unit == "MB" ~ val / 1024,
    unit == "GB" ~ val,
    unit == "TB" ~ val * 1024,
    TRUE ~ NA_real_
  )
}

parse_duration_to_min <- function(x) {
  x <- str_trim(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  h <- suppressWarnings(as.numeric(str_match(x, "(\\d+)\\s*h")[,2]))
  m <- suppressWarnings(as.numeric(str_match(x, "(\\d+)\\s*m")[,2]))
  h[is.na(h)] <- 0
  m[is.na(m)] <- 0
  out <- h * 60 + m
  out[out == 0 & (is.na(x) | x == "")] <- NA_real_
  out
}

normalize_tool <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  x <- str_replace_all(x, "[[:space:]]+", "")
  x <- str_replace_all(x, "_", "")
  x <- str_replace_all(x, "-", "")
  x <- str_replace_all(x, "\\.", "")
  dplyr::case_when(
    str_starts(x, "longtr")  ~ "longtr",
    str_starts(x, "atarva")  ~ "atarva",
    str_starts(x, "medaka")  ~ "medaka",
    str_starts(x, "strkit")  ~ "strkit",
    str_starts(x, "vamos")   ~ "vamos",
    str_starts(x, "strdust") ~ "strdust",
    str_starts(x, "straglr") ~ "straglr",
    TRUE ~ x
  )
}

# -------------------------
# Read table
# -------------------------
df_raw <- if (sep_is_tab) readr::read_tsv(infile, show_col_types = FALSE) else readr::read_csv(infile, show_col_types = FALSE)

df <- df_raw %>%
  rename_with(tolower) %>%
  mutate(tool = normalize_tool(process)) %>%
  filter(tool %in% tools_alpha)

stopifnot(mem_col %in% names(df))
stopifnot("realtime" %in% names(df))

df2 <- df %>%
  mutate(
    mem_gb = parse_mem_to_gb(.data[[mem_col]]),
    rt_min = parse_duration_to_min(realtime)
  )

df_left  <- df2 %>% filter(tool %in% tools_left) %>% mutate(tool = factor(tool, levels = tools_left))
df_strag <- df2 %>% filter(tool == "straglr")   %>% mutate(tool = factor(tool, levels = "straglr"))

# -------------------------
# Common style
# -------------------------
base_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    plot.title   = element_text(face = "bold", size = 13),
    plot.tag     = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.98)
  )

fill_scale_all <- scale_fill_manual(values = color_map, labels = label_map, drop = FALSE)
x_left  <- scale_x_discrete(labels = label_map[tools_left], drop = TRUE)
x_strag <- scale_x_discrete(labels = label_map["straglr"], drop = TRUE)

legend_bottom <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.title = element_blank()
)

# -------------------------
# A: Memory (peak_rss) (tag ONLY on 6-tool plot)
# -------------------------
A_left <- ggplot(df_left, aes(x = tool, y = mem_gb, fill = tool)) +
  geom_boxplot(width = 0.65, outlier.size = 1.1) +
  fill_scale_all + x_left +
  labs(title = "Physical Memory Usage", y = "Memory (GB)", tag = "A") +
  base_theme

A_strag <- ggplot(df_strag, aes(x = tool, y = mem_gb, fill = tool)) +
  geom_boxplot(width = 0.65, outlier.size = 1.1) +
  fill_scale_all + x_strag +
  labs(title = NULL, y = NULL) +
  base_theme +
  theme(plot.tag = element_blank())

figA <- A_left + A_strag + plot_layout(widths = c(6, 1))

# -------------------------
# B: Runtime (tag ONLY on 6-tool plot)
# -------------------------
B_left <- ggplot(df_left, aes(x = tool, y = rt_min, fill = tool)) +
  geom_boxplot(width = 0.65, outlier.size = 1.1) +
  fill_scale_all + x_left +
  labs(title = "Task execution real-time", y = "Execution time (minutes)", tag = "B") +
  base_theme

B_strag <- ggplot(df_strag, aes(x = tool, y = rt_min, fill = tool)) +
  geom_boxplot(width = 0.65, outlier.size = 1.1) +
  fill_scale_all + x_strag +
  labs(title = NULL, y = NULL) +
  base_theme +
  theme(plot.tag = element_blank())

figB <- B_left + B_strag + plot_layout(widths = c(6, 1))

# -------------------------
# STACKED layout: A above B + one shared legend
# -------------------------
final_fig <- (figA / figB) +   # <<< CHANGED (was figA | figB)
  plot_layout(guides = "collect") &
  legend_bottom &
  guides(fill = guide_legend(nrow = 1))

final_fig

# -------------------------
# Save (adjust height because stacked)
# -------------------------
ggsave("HPRC_resource_usage_stacked_AB_peakRSS.png", final_fig, width = 16, height = 7, dpi = 300)
ggsave("HPRC_resource_usage_stacked_AB_peakRSS.pdf", final_fig, width = 16, height = 7)

message("Done. Saved: HPRC_resource_usage_stacked_AB_peakRSS.png and .pdf")