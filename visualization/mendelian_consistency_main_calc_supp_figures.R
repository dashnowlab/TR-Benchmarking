# ============================================================
# MERGED SCRIPT (UPDATED WITH YOUR TOOL PALETTE + LABELS)
#
# Tool palette (your request):
# tableau_colors  <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#FF9DA7")
# tools  <- c("longtr", "atarva", "medaka", "strkit", "vamos", "strdust", "straglr")
# labels <- c("LongTR", "ATaRVa", "Medaka\nTandem", "STRkit", "vamos", "STRdust", "Straglr")
#
# IMPORTANT:
# - Colors are mapped by tool KEYS (not labels) -> robust
# - Display labels are used in legends + axes where relevant
# ============================================================

library(dplyr)
library(stringr)
library(ggpubr)
library(scales)
library(ggplot2)
library(readr)
library(tidyr)
library(forcats)

`%||%` <- function(a, b) if (!is.null(a)) a else b

out_dir <- "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/mendelian_consistency/"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# TABLEAU 10 PALETTE (STATUS COLORS + HEATMAP GRADIENT)
# ============================================================
tableau10 <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
)

# 4-class categorical palette used for Mendelian status
cols <- c(
  "MATCH"     = tableau10[1],  # blue
  "ONE-OFF"   = tableau10[2],  # orange
  "MLEN-OFF"  = tableau10[4],  # teal
  "MISMATCH"  = tableau10[3]   # red
)

# gradient palette for continuous heatmaps
tableau_grad <- c(
  tableau10[1], tableau10[4], tableau10[5],
  tableau10[6], tableau10[2], tableau10[3]
)
tool_levels
# ============================================================
# YOUR TOOL COLORS + DISPLAY LABELS (APPLIED GLOBALLY)
# ============================================================
tableau_colors  <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#FF9DA7")
tools  <- c("longtr", "atarva", "medaka", "strkit", "vamos", "strdust", "straglr")
labels <- c("LongTR", "ATaRVa", "Medaka\nTandem", "STRkit", "vamos", "STRdust", "Straglr")

tool_levels <- tools
tool_labels <- setNames(labels, tools)
tool_cols   <- setNames(tableau_colors, tools)   # <-- map colors by tool KEYS

# ============================================================
# JOBS
# ============================================================
jobs <- list(
  atarva = list(
    dists      = file.path(out_dir, "atarva-dists.txt"),
    catalog    = "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/benchmark-catalog-v2.atarva.bed.gz",
    x5_col     = "X5",
    pos_offset = 1,
    merge_by   = "pos",
    catalog_id_col = "X5"
  ),
  longtr = list(
    dists      = file.path(out_dir, "longtr-dists.txt"),
    catalog    = "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/benchmark-catalog-v2.longtr.bed",
    x5_nchar_col     = "X4",
    pos_offset = 0,
    merge_by   = "id",
    catalog_id_col = "X5"
  ),
  medaka = list(
    dists         = file.path(out_dir, "medaka-dists.txt"),
    catalog       = "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/benchmark-catalog-v2.medaka.bed",
    x5_nchar_col  = "X4",
    pos_offset    = 1,
    merge_by      = "pos",
    catalog_id_col = "X5"
  ),
  straglr = list(
    dists         = file.path(out_dir, "straglr-dists.txt"),
    catalog       = "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/benchmark-catalog-v2.strglr.bed",
    x5_nchar_col  = "X4",
    pos_offset    = 0,
    merge_by      = "pos",
    catalog_id_col = "X4"
  ),
  strdust = list(
    dists         = file.path(out_dir, "strdust-dists.txt"),
    catalog       = "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/benchmark-catalog-v2.strdust.bed",
    x5_nchar_col  = "X4",
    pos_offset    = 1,
    merge_by      = "pos",
    catalog_id_col = "X5"
  ),
  strkit = list(
    dists         = file.path(out_dir, "strkit-dists.txt"),
    catalog       = "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/benchmark-catalog-v2.strkit.bed",
    x5_nchar_col  = "X4",
    pos_offset    = 0,
    merge_by      = "pos",
    catalog_id_col = "X5"
  ),
  vamos = list(
    dists      = file.path(out_dir, "vamos-dists.txt"),
    catalog    = "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/benchmark-catalog-v2.vamos.bed",
    x5_col     = "X7",
    pos_offset = 0,
    merge_by   = "pos",
    catalog_id_col = "X5"
  )
)

# ============================================================
# Helpers
# ============================================================
save_both <- function(p, stem, w, h) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot = p, width = w, height = h, dpi = 300)
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot = p, width = w, height = h, device = "pdf")
}

pick_col <- function(df, base) {
  for (nm in c(base, paste0(base, ".y"), paste0(base, ".x"))) {
    if (nm %in% names(df)) return(nm)
  }
  NA_character_
}

parse_gt_alleles <- function(gt) {
  if (is.na(gt) || gt %in% c("", ".", "./.", ".|.")) return(NA_integer_)
  parts <- str_split(as.character(gt), "[/|]", simplify = TRUE)
  parts <- parts[parts != ""]
  if (length(parts) == 0) return(NA_integer_)
  suppressWarnings(as.integer(parts))
}

calc_kid_dlen_from_seq <- function(ref, alt, gt) {
  if (is.na(ref) || ref == "") return(NA_real_)
  ref_len <- nchar(ref)
  
  alt_vec <- character(0)
  if (!is.na(alt) && alt != "." && alt != "") {
    alt_vec <- str_split(as.character(alt), ",", simplify = FALSE)[[1]]
  }
  
  idx <- parse_gt_alleles(gt)
  if (all(is.na(idx))) return(NA_real_)
  
  allele_seq <- vapply(idx, function(i) {
    if (is.na(i)) return(NA_character_)
    if (i == 0) return(ref)
    if (i <= length(alt_vec)) return(alt_vec[i])
    NA_character_
  }, FUN.VALUE = character(1))
  
  allele_len <- ifelse(is.na(allele_seq), NA_integer_, nchar(allele_seq))
  delta <- abs(allele_len - ref_len)
  out <- suppressWarnings(max(delta, na.rm = TRUE))
  if (is.infinite(out)) NA_real_ else out
}

# ------------------------------------------------------------
# Shared join + motif length extraction
# ------------------------------------------------------------
load_and_join <- function(j) {
  df <- read.delim(j$dists, stringsAsFactors = FALSE)
  
  LEN_COL <- if ("len_dist" %in% names(df)) {
    "len_dist"
  } else if ("kid_dlen" %in% names(df)) {
    "kid_dlen"
  } else {
    NA_character_
  }
  
  if (is.na(LEN_COL)) {
    df <- df %>%
      mutate(
        kid_dlen = mapply(calc_kid_dlen_from_seq, kid_REF, kid_ALT, kid_GT),
        len_dist = kid_dlen
      )
    LEN_COL <- "len_dist"
  }
  
  catalog <- readr::read_tsv(j$catalog, col_names = FALSE, comment = "#", show_col_types = FALSE)
  catalog$X2 <- catalog$X2 + (j$pos_offset %||% 1)
  
  merge_by <- j$merge_by %||% "pos"
  if (merge_by == "pos") {
    df <- df %>% left_join(catalog, by = c("X.chrom" = "X1", "pos" = "X2"))
  } else if (merge_by == "id") {
    if (!("kid_ID" %in% names(df))) stop("merge_by='id' requested but dists has no 'kid_ID'")
    idcol <- j$catalog_id_col %||% "X5"
    if (!(idcol %in% names(catalog))) stop(paste0("merge_by='id' requested but catalog has no column '", idcol, "'"))
    df <- df %>% mutate(kid_ID = as.character(kid_ID))
    catalog[[idcol]] <- as.character(catalog[[idcol]])
    df <- df %>% left_join(catalog, by = setNames(idcol, "kid_ID"))
  } else {
    stop("merge_by must be 'pos' or 'id'")
  }
  
  # motif length extraction (ONE value per row)
  xval <- NA_real_
  if (!is.null(j$x5_col)) {
    nm <- pick_col(df, j$x5_col)
    if (!is.na(nm)) xval <- parse_number(as.character(df[[nm]]))
  } else if (!is.null(j$x5_nchar_col)) {
    nm <- pick_col(df, j$x5_nchar_col)
    if (!is.na(nm)) xval <- nchar(as.character(df[[nm]]))
  }
  
  list(df = df, LEN_COL = LEN_COL, motif_len = xval)
}

# ============================================================
# 1) MENDELIAN CONSISTENCY: per-tool runner (stacked bars)
# ============================================================
run_one_mendel <- function(tool_name, j) {
  
  message("---- MENDEL: ", tool_name, " ----")
  
  lj <- load_and_join(j)
  df <- lj$df
  LEN_COL <- lj$LEN_COL
  x5_val <- lj$motif_len
  
  df_plot <- df %>%
    mutate(
      tool = tool_name,
      kid_GT_raw = as.character(kid_GT),
      kid_GT_grp = case_when(
        kid_GT_raw %in% c("0/0","0|0","0/False") ~ "0/0",
        kid_GT_raw %in% c("0/1","1/0","0|1","1|0") ~ "0/1",
        kid_GT_raw %in% c("1/1","1|1","1/False") ~ "1/1",
        kid_GT_raw %in% c("1/2","2/1","1|2","2|1","2/False") ~ "1/2",
        TRUE ~ NA_character_
      ),
      kid_GT_grp = factor(kid_GT_grp, levels = c("0/0","0/1","1/1","1/2")),
      len_val = parse_number(as.character(.data[[LEN_COL]])),
      x5_val  = x5_val,
      status4 = case_when(
        is.na(len_val) ~ NA_character_,
        len_val == 0 ~ "MATCH",
        (len_val) == 1 ~ "ONE-OFF",
        !is.na(x5_val) & (len_val) <= x5_val ~ "MLEN-OFF",
        TRUE ~ "MISMATCH"
      ),
      status4 = factor(str_squish(status4), levels = names(cols))
    ) %>%
    filter(!is.na(kid_GT_grp), !is.na(status4))
  
  tab_raw <- df_plot %>% count(tool, kid_GT_grp, status4, name = "n")
  
  tab <- tab_raw %>%
    group_by(tool, kid_GT_grp) %>%
    mutate(
      total = sum(n),
      pct = n / total,
      label = ifelse(pct >= 0.04, paste0(percent(pct, accuracy = 1), "\n(n=", n, ")"), "")
    ) %>%
    ungroup()
  
  p <- ggplot(tab, aes(x = kid_GT_grp, y = n, fill = status4)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.8) +
    scale_fill_manual(values = cols, breaks = names(cols), drop = FALSE, name = "Length") +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
    geom_text(aes(label = label),
              position = position_fill(reverse = TRUE, vjust = 0.5),
              size = 3, lineheight = 0.9) +
    geom_text(
      data = tab %>% distinct(tool, kid_GT_grp, total),
      aes(x = kid_GT_grp, y = 1.02, label = paste0("n=", total)),
      inherit.aes = FALSE,
      size = 3
    ) +
    coord_cartesian(ylim = c(0, 1.06), clip = "off") +
    # show pretty tool title
    labs(title = tool_labels[[tool_name]] %||% tool_name, x = "Kid genotype", y = "Percentage of loci") +
    theme_pubr(base_size = 10) +
    theme(
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 8),
      legend.text  = element_text(size = 8),
      plot.margin  = margin(5, 10, 5, 5)
    )
  
  save_both(p, paste0(tool_name, "-mendelian_stacked"), 6.5, 6.0)
  
  list(plot = p, tab = tab, tab_raw = tab_raw)
}

# ------------------------------------------------------------
# Run Mendelian for all tools
# ------------------------------------------------------------
results <- list()
for (tool_name in names(jobs)) {
  results[[tool_name]] <- run_one_mendel(tool_name, jobs[[tool_name]])
}

# ------------------------------------------------------------
# Tables (all tools)
# ------------------------------------------------------------
tab_all <- bind_rows(lapply(results, `[[`, "tab")) %>%
  arrange(tool, kid_GT_grp, status4)

tab_all_raw <- bind_rows(lapply(results, `[[`, "tab_raw")) %>%
  arrange(tool, kid_GT_grp, status4)

write_tsv(tab_all, file.path(fig_dir, "ALL-tools-mendelian_table.tsv"))
write_csv(tab_all, file.path(fig_dir, "ALL-tools-mendelian_table.csv"))
write_tsv(tab_all_raw, file.path(fig_dir, "ALL-tools-mendelian_table_RAW.tsv"))
write_csv(tab_all_raw, file.path(fig_dir, "ALL-tools-mendelian_table_RAW.csv"))

tab_all_wide <- tab_all %>%
  select(tool, kid_GT_grp, status4, n, pct, total) %>%
  pivot_wider(names_from = status4, values_from = c(n, pct), values_fill = 0)

write_tsv(tab_all_wide, file.path(fig_dir, "ALL-tools-mendelian_table_WIDE.tsv"))

# ------------------------------------------------------------
# Big multi-panel stacked plot
# ------------------------------------------------------------
big <- ggplot(tab_all, aes(x = kid_GT_grp, y = n, fill = status4)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.8) +
  scale_fill_manual(values = cols, breaks = names(cols), drop = FALSE, name = "Length") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  geom_text(
    aes(label = label),
    position = position_fill(reverse = TRUE, vjust = 0.5),
    size = 2.8,
    lineheight = 0.9
  ) +
  geom_text(
    data = tab_all %>% distinct(tool, kid_GT_grp, total),
    aes(x = kid_GT_grp, y = 1.02, label = paste0("n=", total)),
    inherit.aes = FALSE,
    size = 2.8
  ) +
  facet_wrap(~ tool, ncol = 3, scales = "fixed",
             labeller = as_labeller(tool_labels)) +
  coord_cartesian(ylim = c(0, 1.06), clip = "off") +
  labs(x = "Kid genotype", y = "Percentage of loci") +
  theme_pubr(base_size = 10) +
  theme(
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 8),
    legend.text  = element_text(size = 8),
    strip.text   = element_text(size = 10, face = "bold"),
    plot.margin  = margin(8, 10, 8, 8)
  )

save_both(big, "ALL-tools-mendelian_stacked_BIG", 12.5, 8.5)

# ============================================================
# 2) EXTRA MENDELIAN PLOTS (ONE GOOD% heatmap only)
# ============================================================

tab_all2 <- tab_all %>%
  mutate(
    tool       = factor(tool, levels = tool_levels),
    kid_GT_grp = factor(kid_GT_grp, levels = c("0/0","0/1","1/1","1/2")),
    status4    = factor(status4, levels = names(cols)),
    pct        = pct %||% (n / total)
  )

good_heat_tab <- tab_all2 %>%
  group_by(tool, kid_GT_grp) %>%
  summarise(
    total    = max(total, na.rm = TRUE),
    good_n   = sum(n[status4 %in% c("MATCH","ONE-OFF","MLEN-OFF")], na.rm = TRUE),
    good_pct = good_n / total,
    .groups  = "drop"
  ) %>%
  mutate(lbl = percent(good_pct, accuracy = 1))

p_good_heat <- ggplot(good_heat_tab, aes(x = kid_GT_grp, y = tool, fill = good_pct)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = lbl), size = 3, lineheight = 0.9) +
  scale_fill_gradientn(
    colors = tableau_grad,
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    name = "GOOD %"
  ) +
  labs(x = "Kid genotype", y = "Tool") +
  scale_y_discrete(limits = tool_levels, labels = tool_labels[tool_levels]) +
  theme_pubr(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_both(p_good_heat, "ALL-tools-GOOD_HEATMAP", 9.5, 5.8)

# ============================================================
# 3) MOTIF-LENGTH ANALYSIS (TABLES + all-tools plot only)
# ============================================================
run_one_motif <- function(tool_name, j) {
  
  message("---- MOTIF: ", tool_name, " ----")
  
  lj <- load_and_join(j)
  df <- lj$df
  LEN_COL <- if ("len_dist" %in% names(df)) "len_dist" else "kid_dlen"
  xval <- lj$motif_len
  
  df2 <- df %>%
    mutate(
      len_val = parse_number(as.character(.data[[LEN_COL]])),
      motif_len = xval,
      status4 = case_when(
        len_val == 0 ~ "MATCH",
        (len_val) == 1 ~ "ONE-OFF",
        (len_val) <= motif_len ~ "MLEN-OFF",
        TRUE ~ "MISMATCH"
      ),
      good = status4 %in% c("MATCH", "ONE-OFF", "MLEN-OFF")
    ) %>%
    filter(!is.na(motif_len))
  
  df2 <- df2 %>%
    mutate(
      motif_len_int = as.integer(round(motif_len)),
      motif_bin = ifelse(motif_len_int > 10, ">10", as.character(motif_len_int)),
      motif_bin = factor(motif_bin, levels = c(as.character(1:10), ">10"))
    )
  
  tab <- df2 %>%
    group_by(motif_bin) %>%
    summarise(
      tool = tool_name,
      n_total = n(),
      n_good = sum(good),
      pct_good = 100 * n_good / n_total,
      .groups = "drop"
    )
  
  tab
}

motif_tabs <- lapply(names(jobs), \(t) run_one_motif(t, jobs[[t]]))
tab_motif_all <- bind_rows(motif_tabs) %>%
  mutate(tool = factor(tool, levels = tool_levels))

write_tsv(tab_motif_all, file.path(fig_dir, "motif_bins_pct_good_all_tools.tsv"))

p_all <- ggplot(tab_motif_all, aes(motif_bin, pct_good, color = tool, group = tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = tool_cols,
    breaks = tool_levels,
    labels = tool_labels[tool_levels],
    drop = FALSE
  ) +
  labs(
    title = "MATCH + ONE-OFF + MLEN-OFF vs motif length (all tools)",
    x = "Motif length bin",
    y = "Percentage (%)",
    color = "Tool"
  ) +
  theme_pubr(base_size = 12)

save_both(p_all, "alltools-motif_line", 11, 5)

# ============================================================
# 4) GOOD% by kid genotype (ALL TOOLS)
# ============================================================
good_tab <- tab_all %>%
  mutate(
    tool = factor(tool, levels = tool_levels),
    kid_GT_grp = factor(kid_GT_grp, levels = c("0/0","0/1","1/1","1/2"))
  ) %>%
  group_by(tool, kid_GT_grp) %>%
  summarise(
    total = max(total, na.rm = TRUE),
    good_n = sum(n[status4 %in% c("MATCH","ONE-OFF","MLEN-OFF")], na.rm = TRUE),
    good_pct = good_n / total,
    .groups = "drop"
  )

write_tsv(good_tab, file.path(fig_dir, "ALL-tools-GOODpct_by_genotype.tsv"))

p_good_alltools <- ggplot(good_tab, aes(x = kid_GT_grp, y = good_pct, group = tool, color = tool)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(
    values = tool_cols,
    breaks = tool_levels,
    labels = tool_labels[tool_levels],
    name = "Tool",
    drop = FALSE
  ) +
  labs(
    x = "Kid genotype",
    y = "GOOD loci (MATCH + ONE-OFF + MLEN-OFF)",
    title = "GOOD rate by kid genotype (all tools)"
  ) +
  theme_pubr(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

save_both(p_good_alltools, "ALL-tools-GOODpct_LINES_alltools", 12.5, 6)

# ============================================================
# ONE-PANEL: GOOD loci COUNTS by kid genotype, per tool
# ============================================================
good_counts <- tab_all %>%
  mutate(
    tool = factor(tool, levels = tool_levels),
    kid_GT_grp = factor(kid_GT_grp, levels = c("0/0","0/1","1/1","1/2"))
  ) %>%
  group_by(tool, kid_GT_grp) %>%
  summarise(
    good_n = sum(n[status4 %in% c("MATCH","ONE-OFF","MLEN-OFF")], na.rm = TRUE),
    total_n = max(total, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(good_counts, file.path(fig_dir, "ALL-tools-GOODcounts_by_genotype.tsv"))

p_good_counts <- ggplot(good_counts, aes(x = kid_GT_grp, y = good_n, fill = tool)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_text(
    aes(label = total_n),
    position = position_dodge(width = 0.85),
    angle = 45,
    vjust = -0.2,
    hjust = 0,
    size = 3
  ) +
  scale_fill_manual(
    values = tool_cols,
    breaks = tool_levels,
    labels = tool_labels[tool_levels],
    name = "Tool",
    drop = FALSE
  ) +
  labs(
    x = "Kid genotype",
    y = "GOOD loci count (MATCH + ONE-OFF + MLEN-OFF)",
    title = "GOOD loci counts by kid genotype (all tools, one panel)"
  ) +
  theme_pubr(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  coord_cartesian(clip = "off")

save_both(p_good_counts, "ALL-tools-GOODcounts_BARS_onepanel", 12.5, 6)

message("DONE. Outputs saved to: ", fig_dir)