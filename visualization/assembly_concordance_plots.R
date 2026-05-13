library(ggplot2)
library(ggridges)
library(cowplot)
library(dplyr)
library(tidyr)
library(stringr)
library(extrafont)
library(ggforce)
library(colorspace)
library(ggnewscale)
library(ggpubr)
library(showtext)
showtext_auto()

font_import()
loadfonts(device = "mac")
fonts()
setwd('/Users/avvarua/Documents/projects/TR-Benchmarking/scripts/')

# color palette for tools
tableau_colors  <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#FF9DA7")
tools <- c("longtr", "atarva", "medaka", "strkit", "vamos", "strdust", "straglr")
labels <- c("LongTR", "ATaRVa", "Medaka\nTandem", "STRkit", "vamos", "STRdust", "Straglr")
tool_labels <- setNames(labels, tools)
tool_colors <- setNames(tableau_colors, labels)
tool_colors["Medaka Tandem"] <- "#E15759"

# Plot themes ------------------------------------------------------------------

font.theme <- theme_minimal(base_family = "Lato")
fontsize.theme <- theme(axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 16),
                        plot.title = element_text(size = 14),
                        strip.text = element_text(size = 12),
                        legend.title = element_text(size = 14),
                        legend.text  = element_text(size = 12),)
legend.no  <- theme(legend.position = "none")
box.theme  <- theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                    axis.ticks.x = element_line(size = 0.4), axis.ticks.length = unit(2, "pt"))

# ------------------------------------------------------------------------------


# Preparing datasets -----------------------------------------------------------
metadata.df <- read.csv(paste0('./hprc-sample-info.tsv'), sep = '\t')
samples_r10 <- (metadata.df %>% filter(metadata.df$flowcell == 'R10.4.1'))$sample_id
samples_r9 <- (metadata.df %>% filter(metadata.df$flowcell == 'R9.4.1'))$sample_id

variant_order <- c('0', '1', '0/0', '0/1', '1/1', '1/2')
variant_order <- c('0/0', '0/1', '1/1', '1/2')
lrange_order <- c("<50bp","50-100bp","100-200bp","200-500bp","500-1000bp",">1000bp")
lrange_labels <- c("Length <50bp", "Length 50-100bp", "Length 100-200bp",
                   "Length 200-500bp", "Length 500-1000bp","Length >1000bp")
mlen_labels <-  c("Motif length: 1bp", "Motif length: 2bp", "Motif length: 3bp", "Motif length: 4bp",
                  "Motif length: 5bp", "Motif length: 6bp", "Motif length: 7bp", "Motif length: 8bp",
                  "Motif length: 9bp", "Motif length: 10bp", "Motif length: >10bp")


combined.df <- NULL
for (i in seq_along(tools)) {
  tool <- tools[i]
  label <- tool_labels[tool]
  print(paste(tool, label))
  df <- read.csv(paste0('hprc-', tool, '-v2.asmcon.tsv'), sep = '\t')
  df$tool <- label
  print(unique(df$tool))
  
  df <- df %>% mutate(tool_allele_bin = cut(tool_len, breaks = c(0, 50, 100, 200, 500, 1000, Inf),
                                            labels = lrange_order, include.lowest = TRUE, right=FALSE))
  df$tool_allele_bin <- factor(df$tool_allele_bin, levels =  lrange_order)
  df <- df %>% mutate(asm_allele_bin = cut(asm_len, breaks = c(0, 50, 100, 200, 500, 1000, Inf),
                                           labels = lrange_order, include.lowest = TRUE, right=FALSE))
  df$asm_allele_bin <- factor(df$asm_allele_bin, levels =  lrange_order)
  df$motif_length_bin <- ifelse(df$motif_length > 10, ">10", as.character(df$motif_length))
  df$motif_length_bin <- factor(df$motif_length_bin, levels = c(as.character(1:10), ">10"), ordered = TRUE)
  df$abs_len_diff <- abs(df$len_diff)
  df$ref_len_diff <- df$asm_len - df$ref_len
  df$abs_ref_len_diff <- abs(df$ref_len_diff)
  df <- df %>% mutate(asm_gt = case_when( asm_vtype == "REF-HOM" ~ "0/0", asm_vtype == "REF-HET" ~ "0/1",
                                          asm_vtype == "ALT-HOM" ~ "1/1", asm_vtype == "ALT-HET" ~ "1/2",
                                          asm_vtype == "REF-HEM" ~ "0",   asm_vtype == "ALT-HEM" ~ "1",
                                          TRUE ~ NA_character_))
  df <- df %>% mutate(tool_gt = case_when( tool_vtype == "REF-HOM" ~ "0/0", tool_vtype == "REF-HET" ~ "0/1",
                                           tool_vtype == "ALT-HOM" ~ "1/1", tool_vtype == "ALT-HET" ~ "1/2",
                                           tool_vtype == "REF-HEM" ~ "0",   tool_vtype == "ALT-HEM" ~ "1",
                                           TRUE ~ NA_character_))
  
  df <- df %>% mutate(asm_gt2 = case_when( asm_vtype == "REF-HOM" ~ "0/0", asm_vtype == "REF-HET" ~ "0/1",
                                           asm_vtype == "ALT-HOM" ~ "1/1", asm_vtype == "ALT-HET" ~ "1/2",
                                           asm_vtype == "REF-HEM" ~ "0/0",   asm_vtype == "ALT-HEM" ~ "1/1",
                                           TRUE ~ NA_character_))
  df <- df %>% mutate(tool_gt2 = case_when( tool_vtype == "REF-HOM" ~ "0/0", tool_vtype == "REF-HET" ~ "0/1",
                                            tool_vtype == "ALT-HOM" ~ "1/1", tool_vtype == "ALT-HET" ~ "1/2",
                                            tool_vtype == "REF-HEM" ~ "0/0",   tool_vtype == "ALT-HEM" ~ "1/1",
                                            TRUE ~ NA_character_))
  df <- df %>% mutate(asm_gt2  = factor(asm_gt2,  levels = variant_order, ordered = TRUE)) %>%
    mutate(tool_gt2 = factor(tool_gt2, levels = variant_order, ordered = TRUE))
  combined.df <- bind_rows(combined.df, df)
}
combined.df$tool <- factor( combined.df$tool, levels = c("LongTR", "ATaRVa", "Medaka\nTandem", "STRkit", "vamos", "STRdust", "Straglr"))
rm(df);


# Preparing datasets -----------------------------------------------------------

# Making Tables ----------------------------------------------------------------

summary.df <- combined.df %>%
  group_by(sample_id, tool, status) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(sample_id, tool, status, fill = list(count = 0)) %>%
  pivot_wider(names_from = status, values_from = count, values_fill = 0) %>%
  mutate(TOTAL = rowSums(across(-c(sample_id, tool)))) %>%
  mutate("MATCH (%)" = if_else(TOTAL == 0, 0, (MATCH/TOTAL)*100),
         "MATCH & ONE-OFF (%)" = if_else(TOTAL == 0, 0, ((MATCH + `ONE-OFF`)/TOTAL)*100),
         "MATCH & ONE-OFF & MOTIF-OFF (%)" = if_else(TOTAL == 0, 0, ((MATCH + `ONE-OFF` + `MOTIF-OFF`)/TOTAL)*100)) %>%
  select(sample_id, tool, TOTAL, MATCH, `ONE-OFF`, `MOTIF-OFF`, MISMATCH, `MATCH (%)`, `MATCH & ONE-OFF (%)`, `MATCH & ONE-OFF & MOTIF-OFF (%)`)%>%
  mutate(tool = ifelse(tool == "Medaka\nTandem", "Medaka Tandem", labels[tool]))
write.table(summary.df, "../plots/hprc-accuracy.tsv", quote = FALSE, row.names = FALSE, sep='\t')

summary.df <- combined.df %>% group_by(sample_id, tool) %>%
  summarise(count = n(), correct = sum(status %in% c("MATCH", "ONE-OFF"), na.rm = TRUE), .groups = "drop") %>%
  mutate(percent = (correct / count) * 100)
stats.df <- summary.df %>% group_by(tool) %>%
  summarise( n_samples = n(), md_percent = median(percent, na.rm = TRUE), mn_percent = mean(percent, na.rm = TRUE),
             sd_percent     = sd(percent, na.rm = TRUE), .groups = "drop")
write.table(stats.df, "../plots/hprc-overall-accuracy.tsv", sep='\t', row.names = FALSE, quote = FALSE,
            col.names = c("Tool", "Samples", "Median Concordance (%)", "Mean Concordance (%)", "SD Concordance"))

summary.df <- combined.df %>%
  group_by(sample_id, catalog_locus) %>%
  summarise(max_allele_length = max(asm_len, na.rm = TRUE), .groups = "drop") %>%
  select(catalog_locus, max_allele_length) %>%
  group_by(catalog_locus) %>%
  summarise(max_allele_length = max(max_allele_length, na.rm = TRUE), .groups = "drop")
write.table(summary.df, "../plots/hprc-locus-max-allele-length.tsv", sep='\t', row.names = FALSE, quote = FALSE)

# Making Tables ----------------------------------------------------------------

filter_samples <- samples_r10

# Supplementart Figure 1 -------------------------------------------------------

summary.nonref.df <- combined.df %>%
  filter(status != "MISSING" & asm_gt2 != "0/0") %>%
  group_by(sample_id, tool) %>%
  summarise(count = n(), correct = sum(status %in% c("MATCH", "ONE-OFF"), na.rm = TRUE), .groups = "drop") %>%
  mutate(percent = (correct / count) * 100)
summary.nonref.df <- merge(summary.nonref.df, metadata.df, by = "sample_id") %>% mutate(Locus = "Variable")
summary.nonref.df <- summary.nonref.df %>% filter(flowcell == "R10.4.1")

summary.df <- combined.df %>%
  filter(status != "MISSING") %>%
  group_by(sample_id, tool) %>%
  summarise(count = n(), correct = sum(status %in% c("MATCH", "ONE-OFF"), na.rm = TRUE), .groups = "drop") %>%
  mutate(percent = (correct / count) * 100)
summary.df <- merge(summary.df, metadata.df, by = "sample_id") %>% mutate(Locus = "All")
summary.df <- summary.df %>% filter(flowcell == "R10.4.1")

summary.df <- bind_rows(summary.df, summary.nonref.df)

plot.df <- summary.df %>% mutate(tool_group = if_else(tool %in% c("LongTR","ATaRVa","Medaka\nTandem","STRkit","vamos"), 
                                                      "Group 1", "Group 2"))

ggplot(plot.df, aes(x = tool, y = percent, fill = tool, alpha = Locus, color = tool)) +
  geom_boxplot( outlier.shape = 21, outlier.color = NA, position = position_dodge(width = 0.6), width = 0.5) +
  scale_x_discrete(drop = TRUE) +
  scale_fill_manual(values = tool_colors, guide = "none") +
  scale_color_manual(values = tool_colors, guide = "none") +
  scale_alpha_manual(values = c( "Variable" = 0.4, "All" = 0.85), 
                     guide = guide_legend(override.aes = list(fill = "#5e5e5e", color = "black"))) +
  labs(x = "",
       y = "Percentage (%)",
       title = "Concordance in All vs Variable loci",
       alpha = "Locus") +
  font.theme + fontsize.theme +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line(size = 0.4),
        axis.ticks.length = unit(2, "pt"),
        strip.text = element_blank())
ggsave('../plots/Supplementary-Figure-1.jpg', width = 6, height = 4)

# Supplementart Figure 1 -------------------------------------------------------

# Supplementary Figure 2 -------------------------------------------------------

filtered.df <- combined.df %>%
  filter(status != "MISSING") %>%
  filter(sample_id %in% samples_r10) %>%
  arrange(asm_gt2)

tool.summary.df <- filtered.df %>%
  select(sample_id, catalog_locus, tool, tool_gt2) %>% distinct() %>%
  group_by(sample_id, tool, tool_gt2) %>%
  summarise(count = n(), .groups="drop") %>%
  group_by(tool, sample_id) %>% mutate(total = sum(count)) %>% ungroup() %>% mutate(percent = (count/total)*100) %>%
  group_by(tool, tool_gt2) %>%
  summarise(avg_count = mean(count), min_count = min(count), max_count = max(count),
            avg_percent = mean(percent), min_percent = min(percent), max_percent = max(percent),.groups = "drop") %>%
  rename("genotype" = "tool_gt2")

asm.summary.df <- filtered.df %>%
  select(sample_id, catalog_locus, tool, asm_gt2) %>% distinct() %>%
  group_by(sample_id, tool, asm_gt2) %>%
  summarise(count = n(), .groups="drop") %>%
  group_by(sample_id, asm_gt2) %>% summarise(count = max(count), .groups = "drop") %>%
  mutate(tool = "ASM") %>%
  group_by(tool, sample_id) %>% mutate(total = sum(count)) %>% ungroup() %>% mutate(percent = (count/total)*100) %>%
  group_by(tool, asm_gt2) %>%
  summarise(avg_count = mean(count), min_count = min(count), max_count = max(count),
            avg_percent = mean(percent), min_percent = min(percent), max_percent = max(percent),.groups = "drop") %>%
  rename("genotype" = "asm_gt2")

summary.df <- rbind(tool.summary.df, asm.summary.df)
summary.df$tool <- factor(summary.df$tool, levels =  c("ASM", "LongTR", "ATaRVa", "Medaka\nTandem", "STRkit", "vamos", "STRdust", "Straglr"))

tool_colors["ASM"] <- "gray50"
ggplot(summary.df, aes(x = tool, y = avg_count, fill = tool, color = tool)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  geom_errorbar( aes(ymin = min_count, ymax = max_count), width = 0.1, position = position_dodge(width = 0.6)) +
  facet_wrap(~genotype, scales = "free") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
  scale_fill_manual(values = tool_colors, guide = "none") +  # <-- custom colors
  scale_color_manual(values = tool_colors, guide = "none") +  # <-- custom colors
  labs( x = "",
    y = "",
    title = "Number of calls in different genotype class") +
  font.theme + fontsize.theme +
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(size = 0.4), axis.ticks.length = unit(2, "pt"),
        axis.text.x = element_text(size = 7))
ggsave('../plots/Supplementary-Figure-2.png', width = 7, height = 4)

# Supplementary Figure 2 -------------------------------------------------------

# Supplementary Figure 3 -------------------------------------------------------

filtered.df <- combined.df %>%
  filter(status != 'MISSING') %>%
  filter(sample_id %in% filter_samples)

tool.summary.df <- filtered.df %>%
  group_by(sample_id, tool, tool_allele_bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(label = paste0("Length ", tool_allele_bin)) %>%
  rename("allele_bin" = tool_allele_bin)%>%
  mutate(tool = ifelse(tool == "Medaka\nTandem", "Medaka Tandem", labels[tool]))
tool.summary.df$allele_bin <- factor(tool.summary.df$allele_bin, levels =  lrange_order)

asm.summary.df <- filtered.df %>%
  select(sample_id, catalog_locus, tool, asm_allele_bin) %>%
  group_by(sample_id, tool, asm_allele_bin) %>%
  summarise(count = n(), .groups="drop") %>%
  group_by(sample_id, asm_allele_bin) %>% summarise(count = max(count), .groups = "drop") %>%
  ungroup() %>%
  mutate(tool = "ASM") %>% 
  mutate(label = paste0("Length ", asm_allele_bin)) %>%
  rename("allele_bin" = asm_allele_bin)
asm.summary.df$allele_bin <- factor(asm.summary.df$allele_bin, levels =  lrange_order)

summary.df <- bind_rows(tool.summary.df, asm.summary.df)
summary.df$tool <- factor(summary.df$tool, levels =  c("ASM", "LongTR", "ATaRVa", "Medaka Tandem", "STRkit", "vamos", "STRdust", "Straglr"))
tool_colors["ASM"] <- "gray50"
tool_colors["Medaka Tandem"] <- tool_colors["Medaka\nTandem"]
summary.df$label <- factor(summary.df$label, levels =  lrange_labels)

plot.df <- summary.df %>%
  group_by(tool, allele_bin, label) %>%
  summarise(avg_count = mean(count),
            min_count = min(count),
            max_count = max(count),
            .groups="drop")
ggplot(plot.df, aes(x = tool, y = avg_count, fill = tool, color = tool)) +
  geom_point(size = 1, position = position_dodge(width = 0.8)) +
  geom_errorbar( aes(ymin = min_count, ymax = max_count), width = 0.1, position = position_dodge(width = 0.8)) +
  facet_wrap(~label, scales = "free", nrow = 2) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
  scale_fill_manual(values = tool_colors) +  # <-- custom colors
  scale_color_manual(values = tool_colors) +  # <-- custom colors
  labs(x = "",
       y = "Number of alleles",
       title = "Number of allele calls per sample - Allele length") +
  font.theme + fontsize.theme +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  guides(color = guide_legend(nrow=2), fill = guide_legend(nrow=2))
ggsave('../plots/Supplementary-Figure-3.jpg', width = 8, height = 6)

# Supplementary Figure 3 -------------------------------------------------------

# Supplementary Figure 4 -------------------------------------------------------

filtered.df <- combined.df %>%
  filter(status != 'MISSING') %>%
  filter(sample_id %in% filter_samples)

summary.df <- filtered.df %>% 
  group_by(tool, sample_id, asm_allele_bin) %>%
  summarise( n_correct = sum(status %in% c("MATCH", "ONE-OFF"), na.rm = TRUE),
             n_motifoff = sum(status %in% c("MATCH", "ONE-OFF", "MOTIF-OFF"), na.rm = TRUE),
             total = n(), .groups = "drop" ) %>%
  group_by(tool, asm_allele_bin) %>%
  summarise( md_percent = median(n_correct, na.rm = TRUE), mn_percent = mean(n_correct, na.rm = TRUE),
             sd_percent     = sd(n_correct, na.rm = TRUE),
             md_mpercent = median(n_motifoff, na.rm = TRUE), mn_mpercent = mean(n_motifoff, na.rm = TRUE),
             sd_mpercent     = sd(n_motifoff, na.rm = TRUE),
             .groups = "drop" )

plot.df <- summary.df %>%
  filter(tool %in% c("LongTR","ATaRVa","Medaka\nTandem","STRkit","vamos","STRdust", "Straglr")) %>%
  mutate(bin_group = if_else(!asm_allele_bin %in% c("500-1000bp",">1000bp"),"Group 1", "Group 2"))

ggplot(plot.df, aes(x = asm_allele_bin, y = mn_percent, color = tool, group = tool)) +
  geom_line(linewidth = 0.6, linetype = "solid", position = position_dodge(width = 0.2)) +
  geom_point(size = 2, alpha = 0.8, position = pos) +
  scale_x_discrete(drop = TRUE) +
  scale_y_log10() +
  scale_color_manual(values = tool_colors) +
  geom_errorbar( aes(ymin = mn_percent - sd_percent, ymax = mn_percent + sd_percent), position = pos, width = 0.1, alpha = 0.8) +
  labs(x = "",
       y = "Number of alleles",
       title = "Accuaracy by counts (Match & Off by 1bp) vs Allele Length Bin",
       subtitle = "Error bars show ±1 SD",
       color = "Tool") +
  font.theme + fontsize.theme +
  theme(legend.position = c(0.85, 0.85),
        strip.text = element_blank())
ggsave('../plots/Supplementary-Figure-4.jpg', width = 9, height = 5)

# Supplementary Figure 4 -------------------------------------------------------

# Supplementary Figure 5 -------------------------------------------------------

filtered.df <- combined.df %>%
  filter(status != "MISSING") %>%
  filter(ref_len_diff > 0) %>%
  filter(sample_id %in% filter_samples)
exp_lrange_order <- c("<50bp","50-100bp","100-200bp","200-500bp","500-1000bp",">1000bp")
filtered.df <- filtered.df %>% mutate(expansion_bin = cut(ref_len_diff, breaks = c(0, 50, 100, 200, 500, 1000, Inf),
                                                          labels = exp_lrange_order, include.lowest = TRUE, right=FALSE))
filtered.df$expansion_bin <- factor(filtered.df$expansion_bin, levels =  exp_lrange_order)

summary.df1 <- filtered.df %>% 
  mutate(tool = 'ASM', "Pore chemistry" = "R10.4.1") %>%
  select(tool, sample_id, catalog_locus, expansion_bin, haplogroup) %>% distinct() %>%
  group_by(tool, expansion_bin) %>% summarise(n_correct = n(), .groups = "drop")

summary.df2 <- filtered.df %>%
  group_by(tool, expansion_bin) %>%
  summarise( n_correct = sum(status %in% c("MATCH", "ONE-OFF", "MOTIF-OFF"), na.rm = TRUE), .groups = "drop" )

summary.df <- merge(summary.df2, summary.df1, by = "expansion_bin") %>%
  select(tool.x, expansion_bin, n_correct.x, n_correct.y) %>%
  rename("Tool" = tool.x, `Expansion bin` = expansion_bin) %>%
  mutate(sensitivity = (n_correct.x/n_correct.y)*100) %>% arrange(Tool, `Expansion bin`) %>%
  filter(`Expansion bin` != "<50bp")

f35 <- ggplot(summary.df, aes(x = `Expansion bin`, y = sensitivity, color = Tool, group = Tool)) +
  geom_line(linewidth = 0) +
  geom_point(size = 3, position = position_dodge(width = 0.4), alpha = 0.7) +
  scale_color_manual(values = tool_colors) +
  labs(x = "Expansion length bin",
       y = "Percentage concordance",
       color = "Tool",
       title = expression("Percentage Concordance (<= Off by motif) of calling expansions")) +
  font.theme +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        strip.text = element_blank())

ggsave('../plots/Supplementary-Figure-5.jpg', width = 6, height = 3.5)

# Supplementary Figure 5 -------------------------------------------------------

# Supplementary Figure 6 -------------------------------------------------------

filtered.df <- combined.df %>%
  filter(status != 'MISSING') %>%
  filter(sample_id %in% filter_samples)

asm.summary.df <- filtered.df %>%
  arrange(motif_length_bin) %>%
  select(sample_id, tool, catalog_locus, motif_length_bin) %>%
  distinct() %>%
  group_by(sample_id, tool, motif_length_bin) %>%
  summarise(count = n(), .groups="drop") %>%
  group_by(sample_id, motif_length_bin) %>% summarise(count = max(count), .groups = "drop") %>%
  ungroup() %>%
  mutate(tool = "ASM") %>% 
  mutate(label = paste0("Motif length: ", motif_length_bin, "bp"))
asm.summary.df$label <- factor(asm.summary.df$label, levels = mlen_labels)

tool.summary.df <- filtered.df %>%
  arrange(motif_length_bin) %>%
  select(sample_id, tool, catalog_locus, motif_length_bin) %>%
  distinct() %>%
  group_by(sample_id, tool, motif_length_bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(label = paste0("Motif length: ", motif_length_bin, "bp")) %>%
  group_by(tool, motif_length_bin) %>%
  summarise(avg_count = mean(count), min_count = min(count), max_count = max(count), .groups = "drop") %>% 
  mutate(label = paste0("Motif length: ", motif_length_bin, "bp"))

tool.summary.df$label <- factor(tool.summary.df$label, levels = mlen_labels)

ggplot(tool.summary.df, aes(x = tool, y = avg_count, fill = tool, color = tool)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  geom_errorbar( aes(ymin = min_count, ymax = max_count), width = 0.1, position = position_dodge(width = 0.6)) +
  facet_wrap(~label, scales = "free_y", nrow = 3) +
  scale_fill_manual(values = tool_colors) +  # <-- custom colors
  scale_color_manual(values = tool_colors) +  # <-- custom colors
  labs(x = "",
       y = "Number of alleles",
       title = "Number of called loci per sample - Motif length bins") +
  font.theme + fontsize.theme +
  theme(legend.position = "none",
        panel.grid.minor =  element_blank(),
        panel.grid.major.x =  element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12))
ggsave('../plots/Supplementary-Figure-6.jpg', width = 12, height = 7)

# Supplementary Figure 6 -------------------------------------------------------

# Figure 3A --------------------------------------------------------------------
summary.df <- combined.df %>%
  filter(status != "MISSING") %>%
  # filter(sample_id %in% filter_samples) %>%
  group_by(sample_id, tool) %>%
  summarise(count = n(), correct = sum(status %in% c("MATCH", "ONE-OFF"), na.rm = TRUE), .groups = "drop") %>%
  mutate(percent = (correct / count) * 100)
summary.df <- merge(summary.df, metadata.df, by = "sample_id") %>% mutate(Locus = "All")

plot.df <- summary.df %>% mutate(tool_group = if_else(tool %in% c("LongTR","ATaRVa","Medaka\nTandem","STRkit","vamos"), 
                                                      "Group 1", "Group 2"))
f3a <- ggplot(plot.df, aes(x = tool, y = percent, fill = tool, alpha = flowcell, color = tool)) +
  geom_boxplot( outlier.shape = 21, outlier.color = NA, position = position_dodge(width = 0.6), width = 0.5) +
  facet_wrap( ~tool_group, scales = "free", space = "free_x", nrow = 1 ) +
  scale_x_discrete(drop = TRUE) +
  scale_fill_manual(values = tool_colors, guide = "none") + 
  scale_color_manual(values = tool_colors, guide = "none") +
  scale_alpha_manual(values = c( "R9.4.1" = 0.4, "R10.4.1" = 0.85), 
                     guide = guide_legend(override.aes = list(fill = "#5e5e5e", color = "black"))) +
  labs(x = "",
       y = "Percentage (%)",
       alpha = "Pore chemistry",
       title = "Percentage Concordance (Match & Off by 1bp)") +
  font.theme + fontsize.theme +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line(size = 0.4),
        axis.ticks.length = unit(2, "pt"),
        strip.text = element_blank())
# Figure 3A --------------------------------------------------------------------

# Figure 3B --------------------------------------------------------------------
combined.df <- combined.df %>% arrange(asm_gt2)
filtered.df <- combined.df %>%
  filter(status != "MISSING") %>%
  filter(sample_id %in% filter_samples)
summary.df <- filtered.df %>%
  group_by(tool, sample_id, asm_gt2) %>%
  summarise( n_correct = sum(status %in% c("MATCH", "ONE-OFF"), na.rm = TRUE), total     = n(), .groups = "drop" ) %>%
  mutate(percent = (n_correct / total) * 100) %>%
  arrange(asm_gt2)
temp.df <- summary.df %>%
  group_by(tool, asm_gt2) %>%
  summarise(avg_percent = mean(percent), median_percent = median(percent), .groups="drop") %>%
  arrange(tool)
pos <- position_jitter(width = 0.2, seed=10)
pos <- position_dodge(width = 0.5)  # adjust width for spacing between tools

# legend_plot <- 

temp.df$tool <- factor( temp.df$tool, levels = c("LongTR", "ATaRVa", "Medaka\nTandem", "STRkit", "vamos", "STRdust", "Straglr"))
temp.df <- temp.df %>%
  mutate(tool = gsub("Medaka\nTandem", "Medaka Tandem", tool)) %>%
  arrange(tool)
temp.df$tool <- factor( temp.df$tool, levels = c("LongTR", "ATaRVa", "Medaka Tandem", "STRkit", "vamos", "STRdust", "Straglr"))
temp.df <- temp.df %>% arrange(tool)

legend_plot <- ggplot(temp.df,
                      aes(x = asm_gt2,
                          y = avg_percent,
                          color = tool)) +
  geom_point(size = 6) +
  scale_color_manual(values = tool_colors, name = NULL) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),      # ⬅ increase text size
    legend.key.size = unit(0.8, "cm")           # ⬅ increase spacing box
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(size = 5)             # ⬅ increase legend dot size
    )
  )

legend_plot <- ggplot(temp.df, aes(x = asm_gt2, y = avg_percent, color = tool)) +
  geom_point(size = 4) +
  scale_color_manual(values = tool_colors, name = NULL) +
  # theme_minimal() +
  theme(legend.position = "bottom") + font.theme + fontsize.theme +
  guides(color = guide_legend(nrow = 1))
legend_plot
legend <- ggpubr::get_legend(legend_plot)

f3b <- ggplot(summary.df, aes(x = asm_gt2, y = percent, color = tool, fill = tool)) +
  geom_boxplot( outlier.shape = 21, outlier.color = NA, position = position_dodge(width = 0.6), width = 0.5, alpha = 0.7) +
  scale_x_discrete(drop = TRUE) +
  scale_fill_manual(values = tool_colors) +  # <-- custom colors
  scale_color_manual(values = tool_colors) +  # <-- custom colors
  font.theme + fontsize.theme +
  labs( x = "Genotype",
        y = "Percentage (%)",
        title = "Percentage Concordance (Match & Off by 1bp) vs Genotype",
        color = "Tool" ) + 
  guides(color = guide_legend(nrow=1), fill = guide_legend(nrow=1)) +
  theme(legend.position = "none", legend.title = element_blank())
f3b <- f3b + theme(legend.position = "none")

# Figure 3B --------------------------------------------------------------------

# Figure 3C --------------------------------------------------------------------

combined.df <- combined.df %>% arrange(asm_allele_bin)
filtered.df <- combined.df %>%
  filter(status != "MISSING") %>%
  filter(sample_id %in% filter_samples)
summary.df <- filtered.df %>% 
  group_by(tool, sample_id, asm_allele_bin) %>%
  summarise( n_correct = sum(status %in% c("MATCH", "ONE-OFF"), na.rm = TRUE),
             total     = n(), .groups = "drop" ) %>%
  mutate(percent = (n_correct / total) * 100) %>%
  group_by(tool, asm_allele_bin) %>%
  summarise( md_percent = median(percent, na.rm = TRUE),
             mn_percent = mean(percent, na.rm = TRUE),
             sd_percent     = sd(percent, na.rm = TRUE), .groups = "drop" )

pos <- position_dodge(width = 0.2)
plot.df <- summary.df %>%
  filter(tool %in% c("LongTR","ATaRVa","Medaka\nTandem","STRkit","vamos")) %>%
  mutate(bin_group = if_else(!asm_allele_bin %in% c("500-1000bp",">1000bp"),"Group 1", "Group 2"))

p1 <- ggplot(plot.df, aes(x = asm_allele_bin, y = mn_percent, color = tool, group = tool)) +
  geom_line(linewidth = 0.6, linetype = "solid", position = pos) +
  geom_point(size = 2, alpha = 0.8, position = pos) +
  # facet_wrap(~bin_group, scales = "free", space = "free_x") +
  scale_x_discrete(drop = TRUE) +
  scale_color_manual(values = tool_colors) +  # <-- custom colors
  geom_errorbar( aes(ymin = mn_percent - sd_percent, ymax = mn_percent + sd_percent), position = pos, width = 0.1, alpha = 0.8) +
  labs(x = "", y = "Percentage (%)", title = expression("Percentage (<= Off by motif) vs Allele Length Bin"),
       subtitle = "Error bars show ±1 SD", color = "Tool") +
  font.theme + fontsize.theme +
  theme(axis.text.x = element_text(size = 13.6)) +
  theme(strip.text = element_blank(), legend.position = "none")

plot.df <- summary.df %>% filter(!tool %in% c("LongTR","ATaRVa","Medaka\nTandem","STRkit","vamos"))
p2 <- ggplot(plot.df, aes(x = asm_allele_bin, y = mn_percent, color = tool, group = tool)) +
  geom_line(linewidth = 0.6, linetype = "solid", position = pos) +
  geom_point(size = 2, alpha = 0.8, position = pos) +
  # facet_wrap(~bin_group, scales = "free", space = "free_x") +
  scale_x_discrete(drop = TRUE) +
  scale_color_manual(values = tool_colors) +  # <-- custom colors
  geom_errorbar( aes(ymin = mn_percent - sd_percent, ymax = mn_percent + sd_percent), position = pos, width = 0.1, alpha = 0.8) +
  labs(x = "Allele length bin", y = "Percentage (%)", color = "Tool") +
  font.theme + fontsize.theme +
  theme(axis.text.x = element_text(size = 13.6)) +
  theme(strip.text = element_blank(), legend.position = "none")
f3c <- plot_grid(p1, p2, ncol = 1, rel_heights = c(2,1))

# Figure 3C --------------------------------------------------------------------

# Figure 3D --------------------------------------------------------------------

filtered.df <- combined.df %>%
  filter(status != 'MISSING') %>%
  filter(sample_id %in% filter_samples)
summary.df <- filtered.df %>%
  group_by(tool, sample_id, motif_length_bin) %>%
  summarise(correct = sum(status %in% c("MATCH", "ONE-OFF"), na.rm = TRUE), count = n(), .groups = "drop" ) %>%
  mutate(percent = (correct / count) * 100) %>%
  group_by(tool, motif_length_bin) %>%
  summarise( median_percent = median(percent, na.rm = TRUE), mean_percent = mean(percent, na.rm = TRUE),
             sd_percent = sd(percent, na.rm = TRUE), .groups = "drop") %>% 
  mutate(label = paste0("Motif length: ", motif_length_bin, "bp"))

plot.df <- summary.df %>% filter(tool %in% c("LongTR","ATaRVa","Medaka\nTandem","STRkit","vamos")) %>%
  mutate(bin_group = if_else(motif_length_bin <= 2,"Group 1", "Group 2"))
p1 <- ggplot(plot.df, aes(x = motif_length_bin, y = mean_percent, color = tool, group = tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1, alpha = 0.8) +
  facet_wrap(~bin_group, scales = "free", nrow = 1) +
  geom_errorbar( aes(ymin = mean_percent - sd_percent, ymax = mean_percent + sd_percent), width = 0.1, alpha = 0.8) +
  scale_color_manual(values = tool_colors) +
  labs(x = "", y = "Percentage (%)", title = "Percentage (Match & Off by 1bp) vs Motif length",
       subtitle = "Error bars show ±1 SD", color = "Tool") +
  font.theme + fontsize.theme +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                     strip.text = element_blank(), legend.position = "none")

plot.df <- summary.df %>% filter(!tool %in% c("LongTR","ATaRVa","Medaka\nTandem","STRkit","vamos")) %>%
  mutate(bin_group = if_else(motif_length_bin <= 3,"Group 1", "Group 2"))
p2 <- ggplot(plot.df, aes(x = motif_length_bin, y = mean_percent, color = tool, group = tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1, alpha = 0.8) +
  # facet_wrap(~bin_group, scales = "free", space = "free_x") +
  geom_errorbar( aes(ymin = mean_percent - sd_percent, ymax = mean_percent + sd_percent), width = 0.1, alpha = 0.8) +
  scale_color_manual(values = tool_colors) +
  labs(x = "Motif length (bp)", y = "Percentage (%)", color = "Tool") +
  font.theme + fontsize.theme +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                     strip.text = element_blank(), legend.position = "none")
f3d <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1.5,1))

# Figure 3D --------------------------------------------------------------------

# Figure 3E ------------------------------------------------------

filtered.df <- combined.df %>%
  filter(status != "MISSING") %>%
  filter(sample_id %in% filter_samples)

summary.df1 <- filtered.df %>%
  filter(tool != "Straglr") %>%
  select(sample_id, tool, status) %>%
  group_by(sample_id, tool, status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample_id, tool) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(percent = (count/total)*100 ) %>%
  filter(status == "MATCH") %>%
  mutate(match = "Length")
  
summary.df2 <- filtered.df %>%
  filter(tool != "Straglr") %>%
  select(sample_id, tool, tool_lev) %>%
  group_by(sample_id, tool, tool_lev) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample_id, tool) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(percent = (count/total)*100 ) %>%
  filter(tool_lev == 0) %>%
  mutate(match = "Lev Dist")

summary.df = bind_rows(summary.df1, summary.df2)

stats.df <- summary.df %>%
  group_by(tool, match) %>%
  summarise(average = mean(percent), .groups = "drop") %>%
  pivot_wider(names_from = match, values_from = average) %>%
  mutate(difference = Length - `Lev Dist`)

plot.df <- summary.df %>% mutate(tool_group = if_else(tool %in% c("LongTR","ATaRVa","Medaka\nTandem","STRkit","vamos"), 
                                                      "Group 1", "Group 2"))

f3e <- ggplot(plot.df, aes(x = tool, y = percent, fill = tool, alpha = match, color = tool)) +
  geom_boxplot( outlier.shape = 21, outlier.color = NA, position = position_dodge(width = 0.6), width = 0.5) +
  facet_wrap(~tool_group, scales = "free", space = "free_x") +
  scale_x_discrete(drop = TRUE) +
  scale_fill_manual(values = tool_colors, guide = "none") +
  scale_color_manual(values = tool_colors, guide = "none") +
  scale_alpha_manual(values = c( "Lev Dist" = 0.4, "Length" = 0.85), 
                     guide = guide_legend(override.aes = list(fill = "#5e5e5e", color = "black"))) +
  labs(x = "",
       y = "Percentage (%)",
       title = "Concordance based on Length and Sequence",
       alpha = "Concordance type") +
  font.theme + fontsize.theme +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line(size = 0.4),
        axis.ticks.length = unit(2, "pt"),
        strip.text = element_blank())
# Figure 3E --------------------------------------------------------------------

# Figure 3F --------------------------------------------------------------------

filtered.df <- combined.df %>%
  filter(status != 'MISSING') %>%
  filter(sample_id %in% filter_samples) %>%
  filter(tool != 'Straglr') %>%
  mutate(lev_dist = tool_lev - abs_len_diff)

totals.df <- filtered.df %>%
  group_by(sample_id, tool) %>%
  summarise(total = n(), .groups = "drop")

dist.df <- filtered.df %>%
  filter(status == "MATCH")

limit <- 5
limit_label <- paste0(">", limit)
dist.df$lev_bin <- ifelse(dist.df$lev_dist > limit, limit_label, as.character(dist.df$lev_dist))
dist.df$lev_bin <- factor(dist.df$lev_bin, levels = c(as.character(0:10), limit_label), ordered = TRUE)

dist.df <- dist.df %>%
  group_by(tool, sample_id, lev_bin) %>%
  summarise(count = n(), .groups = "drop") #%>%
dist.df <- merge(dist.df, totals.df, by=c("tool", "sample_id"))
dist.df <- dist.df %>%
  mutate(percent = (count/total)*100) %>%
  group_by(tool, lev_bin) %>%
  summarise(average = mean(percent), .groups = "drop")

stats.df <- dist.df %>%
  pivot_wider(names_from = lev_bin, values_from = average) %>%
  mutate(sum = `1` + `2`)

plot.df <- dist.df %>%
  filter(lev_bin > 0) %>%
  mutate(group = ifelse(tool == "STRdust", "Group2", "Group1"))
f3f <- ggplot(plot.df, aes(x = lev_bin, y = average, fill = tool, color = tool)) +
  geom_col(alpha = 0.7, position = "identity") +
  scale_fill_manual(values = tool_colors) +
  scale_color_manual(values = tool_colors) +
  facet_wrap(~ tool, ncol = 3) +
  theme(legend.position = "none") +
  labs( y = "Percentage(%)", x = "Levenshtein distance", title = "Concordance drop at varying levenshtein distances" ) +
  font.theme + fontsize.theme +
  theme( legend.position = "none",
    strip.text = element_text(size = 10, lineheight = 1.1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.ticks.x = element_line(size = 0.4),
    axis.ticks.length = unit(2, "pt"))

# Figure 3F --------------------------------------------------------------------

row_1 = plot_grid(f3a, f3b, labels = c ("A", "B"), ncol = 2)
row_2 <- plot_grid(NULL, legend, NULL, ncol = 3,
                   rel_widths = c(0.75, 1, 1)   # middle column controls legend width
)
row_3 = plot_grid(f3c, f3d, labels = c ("C", "D"), ncol = 2)
row_4 = plot_grid(f3e, f3f, labels = c ("E", "F"), ncol = 2)
final_plot <- plot_grid(
  row_1,
  row_2,
  row_3,
  row_4,
  ncol = 1,
  rel_heights = c(1, 0.2, 1, 1)  # make legend row shorter
)
final_plot
ggsave('../plots/Figure-3.pdf', width = 14, height = 16)

row_1 = plot_grid(f3b, f3c, labels = c ("A", "B"), ncol = 2)
row_2 <- plot_grid(NULL, legend, NULL, ncol = 3,
  rel_widths = c(0.75, 1, 1)   # middle column controls legend width
)
row_3 = plot_grid(f3d, f3e, labels = c ("C", "D"), ncol = 2)
row_4 = plot_grid(f3f, f35, labels = c ("E", "F"), ncol = 2)
final_plot <- plot_grid(
  row_1,
  row_2,
  row_3,
  row_4,
  ncol = 1,
  rel_heights = c(1, 0.2, 1, 1)  # make legend row shorter
)
final_plot
ggsave('../plots/Figure-3R9-fontfix.jpg', p = final_pot, width = 14, height = 16)
