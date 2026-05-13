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


# Figure 5 ---------------------------------------------------------------------
tools_heatmap  <- c("ASM", "LongTR", "Atarva", "Medaka.Tandem", "Vamos", "STRkit", "STRDust", "Straglr")
labels_heatmap <- c("ASM", "LongTR", "ATaRVa", "Medaka\nTandem", "vamos", "STRkit", "STRdust", "Straglr")

matrix.df <- NULL
for (i in seq_along(samples_r10)) {
  sample <- samples_r10[i]
  # if (sample != 'HG00290') { next; }
  df <- read.csv(paste0('../data/compare-tools/', sample, '-lev-comp.tsv.gz'), sep = '\t', comment.char = '#') %>%
    mutate(LENGTH = END-START)
  matrix.df <- bind_rows(matrix.df, df)
}


matrix.df <- matrix.df %>% mutate(across(where(is.numeric), abs))
col.means <- matrix.df %>% summarise(across(6:ncol(matrix.df)-1, \(x)mean(x, na.rm = TRUE)))
out <- col.means %>% pivot_longer(everything(), names_to = "comparison", values_to = "value") %>%
  rowwise() %>%
  mutate( Tool_1 = tools_heatmap[sapply(tools_heatmap, function(t) startsWith(comparison, paste0(t, ".")))][1],
          Tool_2 = sub(paste0("^", Tool_1, "\\."), "", comparison) ) %>%
  ungroup() %>% select(Tool_1, Tool_2, value)
heat.lev.df <- bind_rows(out, out %>% rename(Tool_1 = Tool_2, Tool_2 = Tool_1))

matrix.df <- NULL
for (i in seq_along(samples_r10)) {
  sample <- samples_r10[i]
  df <- read.csv(paste0('../data/compare-tools/', sample, '-len-comp.tsv.gz'), sep = '\t', comment.char = '#') %>%
    mutate(LENGTH = END-START)
  matrix.df <- bind_rows(matrix.df, df)
}

matrix.df <- matrix.df %>% mutate(across(where(is.numeric), abs))
col.means <- matrix.df %>% summarise(across(6:ncol(matrix.df)-1, \(x)mean(x, na.rm = TRUE)))
out <- col.means %>% pivot_longer(everything(), names_to = "comparison", values_to = "value") %>%
  rowwise() %>%
  mutate( Tool_1 = tools_heatmap[sapply(tools_heatmap, function(t) startsWith(comparison, paste0(t, ".")))][1],
          Tool_2 = sub(paste0("^", Tool_1, "\\."), "", comparison) ) %>%
  ungroup() %>% select(Tool_1, Tool_2, value)
heat.len.df <- bind_rows(out, out %>% rename(Tool_1 = Tool_2, Tool_2 = Tool_1))

heat.lev.df <- heat.lev.df %>% mutate( Tool_1 = factor(Tool_1, levels = tools_heatmap), Tool_2 = factor(Tool_2, levels = tools_heatmap) )
heat.len.df <- heat.len.df %>% mutate( Tool_1 = factor(Tool_1, levels = tools_heatmap), Tool_2 = factor(Tool_2, levels = tools_heatmap) )

diag_df <- data.frame( Tool_1 = tools_heatmap, Tool_2 = tools_heatmap, value = 0 )
heat.lev.upper <- heat.lev.df %>% filter(as.numeric(Tool_1) <= as.numeric(Tool_2)) %>%
  bind_rows(diag_df)  %>% mutate( Tool_1 = factor(Tool_1, levels = tools_heatmap), Tool_2 = factor(Tool_2, levels = tools_heatmap) ) %>% arrange(Tool_1, Tool_2)
heat.len.lower <- heat.len.df %>% filter(as.numeric(Tool_1) > as.numeric(Tool_2))%>%
  bind_rows(diag_df) %>% mutate( Tool_1 = factor(Tool_1, levels = tools_heatmap), Tool_2 = factor(Tool_2, levels = tools_heatmap) ) %>% arrange(Tool_1, Tool_2)


ggplot() +
  # upper triangle (df1)
  geom_tile( data = heat.lev.upper, aes(x=Tool_1, y=Tool_2, fill = value), color = "white" ) +
  geom_text(data=heat.lev.upper, aes(x=Tool_1, y=Tool_2, label = round(value, 2)), size = 12) +
  scale_fill_gradient( low = "white", high = "#D73027", name = "Avg Lev. Dist" ) +
  
  ggnewscale::new_scale_fill() +

  # lower triangle (df2)
  geom_tile( data = heat.len.lower, aes(x=Tool_1, y=Tool_2, fill = value), color = "white" ) +
  geom_text(data=heat.len.lower, aes(x=Tool_1, y=Tool_2, label = round(value, 2)), size = 12) +
  scale_fill_gradient( low = "white", high = "#4575B4", name = "Avg Length Diff" ) +
  scale_x_discrete(labels=labels_heatmap) +
  scale_y_discrete(labels=labels_heatmap) +
  coord_equal() + font.theme + fontsize.theme +
  theme(axis.text.x = element_text(size = 36),
        axis.text.y = element_text(size = 39),
        axis.title.x = element_text(size = 48),
        axis.title.y = element_text(size = 48),
        plot.title = element_text(size = 48),
        legend.title = element_text(size = 42),
        legend.text  = element_text(size = 42),) +
  labs(x = "", y = "", title =  "Comparison of Tools") +
  theme( axis.text.x = element_text(hjust = 0.5), panel.grid.major = element_blank())

ggsave('../plots/Figure-5.png', width = 4, height = 3.5)
# ggsave('../plots/Supplementary-Figure-9.jpg', width = 8, height = 7)

# Figure 5 --------------------------------------------------------------------- 

# Supplementary Figure 9 -------------------------------------------------------

tools_heatmap  <- c("ASM", "LongTR", "Atarva", "Medaka.Tandem", "Vamos", "STRkit", "STRdust", "Straglr")
labels_heatmap <- c("ASM", "LongTR", "ATaRVa", "Medaka\nTandem", "vamos", "STRkit", "STRdust", "Straglr")

regions.df <- read.csv('../scripts/HG00290-mismatch.txt', sep = '\t', comment.char = '#')

matrix.df <- NULL
for (i in seq_along(samples_r10)) {
  sample <- samples_r10[i]
  if (sample != 'HG00290') { next; }
  df <- read.csv(paste0('../data/compare-tools/', sample, '-lev-comp.tsv.gz'), sep = '\t', comment.char = '#') %>%
    mutate(LENGTH = END-START)
  matrix.df <- bind_rows(matrix.df, df)
}

matrix.df <- matrix.df %>% semi_join(regions.df, by=c("CHROM", "START", "END"))
matrix.df <- matrix.df %>% mutate(across(where(is.numeric), abs))
col.means <- matrix.df %>% summarise(across(6:ncol(matrix.df)-1, \(x)mean(x, na.rm = TRUE)))
out <- col.means %>% pivot_longer(everything(), names_to = "comparison", values_to = "value") %>%
  rowwise() %>%
  mutate( Tool_1 = tools_heatmap[sapply(tools_heatmap, function(t) startsWith(comparison, paste0(t, ".")))][1],
          Tool_2 = sub(paste0("^", Tool_1, "\\."), "", comparison) ) %>%
  ungroup() %>% select(Tool_1, Tool_2, value)
heat.lev.df <- bind_rows(out, out %>% rename(Tool_1 = Tool_2, Tool_2 = Tool_1))

matrix.df <- NULL
for (i in seq_along(samples_r10)) {
  sample <- samples_r10[i]
  if (sample != 'HG00290') { next; }
  df <- read.csv(paste0('../data/compare-tools/', sample, '-len-comp.tsv.gz'), sep = '\t', comment.char = '#') %>%
    mutate(LENGTH = END-START)
  matrix.df <- bind_rows(matrix.df, df)
}

matrix.df <- matrix.df %>% semi_join(regions.df, by=c("CHROM", "START", "END"))
matrix.df <- matrix.df %>% mutate(across(where(is.numeric), abs))
col.means <- matrix.df %>% summarise(across(6:ncol(matrix.df)-1, \(x)mean(x, na.rm = TRUE)))
out <- col.means %>% pivot_longer(everything(), names_to = "comparison", values_to = "value") %>%
  rowwise() %>%
  mutate( Tool_1 = tools_heatmap[sapply(tools_heatmap, function(t) startsWith(comparison, paste0(t, ".")))][1],
          Tool_2 = sub(paste0("^", Tool_1, "\\."), "", comparison) ) %>%
  ungroup() %>% select(Tool_1, Tool_2, value)
heat.len.df <- bind_rows(out, out %>% rename(Tool_1 = Tool_2, Tool_2 = Tool_1))

heat.lev.df <- heat.lev.df %>% mutate( Tool_1 = factor(Tool_1, levels = tools_heatmap), Tool_2 = factor(Tool_2, levels = tools_heatmap) )
heat.len.df <- heat.len.df %>% mutate( Tool_1 = factor(Tool_1, levels = tools_heatmap), Tool_2 = factor(Tool_2, levels = tools_heatmap) )

diag_df <- data.frame( Tool_1 = tools_heatmap, Tool_2 = tools_heatmap, value = 0 )
heat.lev.upper <- heat.lev.df %>% filter(as.numeric(Tool_1) <= as.numeric(Tool_2)) %>%
  bind_rows(diag_df)  %>% mutate( Tool_1 = factor(Tool_1, levels = tools_heatmap), Tool_2 = factor(Tool_2, levels = tools_heatmap) ) %>% arrange(Tool_1, Tool_2)
heat.len.lower <- heat.len.df %>% filter(as.numeric(Tool_1) > as.numeric(Tool_2))%>%
  bind_rows(diag_df) %>% mutate( Tool_1 = factor(Tool_1, levels = tools_heatmap), Tool_2 = factor(Tool_2, levels = tools_heatmap) ) %>% arrange(Tool_1, Tool_2)


ggplot() +
  # upper triangle (df1)
  geom_tile( data = heat.lev.upper, aes(x=Tool_1, y=Tool_2, fill = value), color = "white" ) +
  geom_text(data=heat.lev.upper, aes(x=Tool_1, y=Tool_2, label = round(value, 2)), size = 4) +
  scale_fill_gradient( low = "white", high = "#D73027", name = "Avg Lev. Dist" ) +
  
  ggnewscale::new_scale_fill() +
  
  # lower triangle (df2)
  geom_tile( data = heat.len.lower, aes(x=Tool_1, y=Tool_2, fill = value), color = "white" ) +
  geom_text(data=heat.len.lower, aes(x=Tool_1, y=Tool_2, label = round(value, 2)), size = 4) +
  scale_fill_gradient( low = "white", high = "#4575B4", name = "Avg Length Diff" ) +
  scale_x_discrete(labels=labels_heatmap) +
  scale_y_discrete(labels=labels_heatmap) +
  coord_equal() + font.theme + fontsize.theme +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 13),) +
  labs(x = "", y = "", title =  "Comparison of Tools") +
  theme( axis.text.x = element_text(hjust = 0.5), panel.grid.major = element_blank())

ggsave('../plots/Supplementary-Figure-9.png', width = 8, height = 7)

# Supplementary Figure 9 -------------------------------------------------------