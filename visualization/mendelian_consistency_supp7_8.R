suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

out_dir <- "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/mendelian_consistency/figures/"
mendel_tsv <- file.path(out_dir, "ALL-tools-mendelian_table_WIDE.tsv")

tool_labels <- c(
  atarva  = "ATaRVa",
  longtr  = "LongTR",
  medaka  = "MedakaTandem",
  straglr = "Straglr",
  strdust = "STRdust",
  strkit  = "STRkit",
  vamos   = "Vamos"
)

pct_fmt <- function(x, digits = 2) sprintf(paste0("%.", digits, "f%%"), x)

mendel <- read_tsv(mendel_tsv, show_col_types = FALSE) %>%
  mutate(
    tool = tolower(tool),
    Tool = tool_labels[tool],
    kid_GT_grp = factor(kid_GT_grp, levels = c("0/0","0/1","1/1","1/2"))
  )

supp7 <- mendel %>%
  transmute(
    Tool,
    kid_GT_grp,
    Match      = 100 * n_MATCH / total,
    `1bp off`  = 100 * `n_ONE-OFF` / total,
    `mlen off` = 100 * `n_MLEN-OFF` / total,
    Mismatch   = 100 * n_MISMATCH / total,
    Total      = 100 * (n_MATCH + `n_ONE-OFF` + `n_MLEN-OFF`) / total
  ) %>%
  pivot_longer(
    cols = c(Match, `1bp off`, `mlen off`, Mismatch, Total),
    names_to = "Metric",
    values_to = "pct"
  ) %>%
  mutate(pct = pct_fmt(pct)) %>%
  pivot_wider(names_from = kid_GT_grp, values_from = pct) %>%
  arrange(Tool, factor(Metric,
                       levels = c("Match","1bp off","mlen off","Mismatch","Total")))

write_tsv(supp7, file.path(out_dir, "Supplementary_Table_7_HG002_trio.tsv"))

# Optional Excel export
wb <- createWorkbook()
addWorksheet(wb, "SuppTable7_GenotypeClass")
writeData(wb, "SuppTable7_GenotypeClass", supp7)
saveWorkbook(wb,
             file.path(out_dir, "Supplementary_Table_7_HG002_trio.xlsx"),
             overwrite = TRUE)



# ============================================================
# Supplementary Table 8 | Mendelian Consistency by Motif Length — HG002 Trio
# (GOOD = MATCH + ±1bp + <1 motif unit)
# Input: motif_bins_pct_good_all_tools.tsv  (from your script)
# Output: Supplementary_Table_8_HG002_trio.tsv + .xlsx
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
})

out_dir  <- "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/mendelian_consistency/figures/"
motif_tsv <- file.path(out_dir, "motif_bins_pct_good_all_tools.tsv")

tool_labels <- c(
  atarva  = "ATaRVa",
  longtr  = "LongTR",
  medaka  = "MedakaTandem",
  straglr = "Straglr",
  strdust = "STRdust",
  strkit  = "STRkit",
  vamos   = "Vamos"
)

# tool order for the table (edit if you want a different order)
tool_order <- c("atarva","longtr","medaka","straglr","strdust","strkit","vamos")

pct_fmt <- function(x, digits = 2) sprintf(paste0("%.", digits, "f%%"), x)

motif <- read_tsv(motif_tsv, show_col_types = FALSE)

# ---- Handle BOTH possible input shapes ----
# A) long format: tool, motif_bin, pct_good
# B) wide format: tool + columns for motif bins
if (all(c("tool","motif_bin","pct_good") %in% names(motif))) {
  
  supp8 <- motif %>%
    mutate(
      tool = tolower(tool),
      Tool = unname(tool_labels[tool]),
      tool = factor(tool, levels = tool_order),
      motif_bin = as.character(motif_bin),
      # make screenshot-like headers
      motif_col = case_when(
        motif_bin == ">10" ~ "Motiflength:>10bp",
        TRUE ~ paste0("Motiflength:", motif_bin, "bp")
      ),
      motif_col = factor(motif_col,
                         levels = c(paste0("Motiflength:", 1:10, "bp"), "Motiflength:>10bp"))
    ) %>%
    transmute(Tool, tool, motif_col, pct = pct_fmt(pct_good)) %>%
    arrange(tool, motif_col) %>%
    select(-tool) %>%
    pivot_wider(names_from = motif_col, values_from = pct)
  
} else {
  
  # wide format: rename motif columns into Motiflength:*bp
  stopifnot("tool" %in% names(motif))
  
  bin_cols <- setdiff(names(motif), "tool")
  
  supp8 <- motif %>%
    mutate(
      tool = tolower(tool),
      Tool = unname(tool_labels[tool]),
      tool = factor(tool, levels = tool_order)
    ) %>%
    arrange(tool) %>%
    mutate(across(all_of(bin_cols), ~ pct_fmt(as.numeric(.x)))) %>%
    select(Tool, all_of(bin_cols)) %>%
    rename_with(function(x) {
      if (x == "tool" || x == "Tool") return(x)
      if (x == ">10") return("Motiflength:>10bp")
      if (grepl("^\\d+$", x)) return(paste0("Motiflength:", x, "bp"))
      # keep if already Motiflength:*
      x
    }, all_of(bin_cols))
}

# Write outputs
write_tsv(supp8, file.path(out_dir, "Supplementary_Table_8_HG002_trio.tsv"))

wb <- createWorkbook()
addWorksheet(wb, "SuppTable8_MotifLength")
writeData(wb, "SuppTable8_MotifLength", supp8)
saveWorkbook(wb, file.path(out_dir, "Supplementary_Table_8_HG002_trio.xlsx"), overwrite = TRUE)

cat("Wrote:\n",
    file.path(out_dir, "Supplementary_Table_8_HG002_trio.tsv"), "\n",
    file.path(out_dir, "Supplementary_Table_8_HG002_trio.xlsx"), "\n")