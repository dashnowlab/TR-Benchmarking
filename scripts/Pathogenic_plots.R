library(ggplot2)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot()) # Sets the default for subsequent plots

tool_colours = c("#7FB3D5", "#F5A623", "#7DCEA0", "#F1948A", "#B39DDB", "#F7DC6F")

data <- read_tsv("pathogenic_results.tsv")
data$tool_allele_max = pmax(data$allele1_len, data$allele2_len, na.rm = T)
data$molec_allele_max = pmax(data$molec_allele1_len, data$molec_allele2_len, na.rm = T)
data$allele_max_diff = data$tool_allele_max - data$molec_allele_max
data$gene_sample = factor(interaction(data$gene, data$sample, sep = " | ", drop = TRUE))
data$molec_type[data$molec_type == "Unknown"] = NA

## Pathogenic loci "strip" plot

# Positions of shaded rectangles
pathogenic_ranges <- data.frame(
  gene = c("HTT", "RFC1", "ATXN1", 
                  "C9orf72", "DMPK", "FMR1",
                  "FXN", "NOTCH2NLC", "PABPN1"),
  path_min = c(36, 400, 39,
           31, 50, 201,
           56, 66, 12),
  path_max = c(250, 2750, 91,
           4088, 4000, 2000,
           1700, 517, 18),
  motif_size = c(3, 5, 3,
                 6, 3, 3,
                 3, 3, 3)
)
pathogenic_ranges$xmin = pathogenic_ranges$path_min * pathogenic_ranges$motif_size
pathogenic_ranges$xmax = pathogenic_ranges$path_max * pathogenic_ranges$motif_size

gene_max_x <- data %>%
  summarise(
    max_x_val = max(tool_allele_max, molec_allele_max, na.rm = TRUE),
    .by = gene
  )

# Set max point on plot
# max_x_val = max(data$tool_allele_max, data$molec_allele_max, na.rm = T)
pathogenic_ranges = merge(pathogenic_ranges, gene_max_x)
pathogenic_ranges$xmax = pmin(pathogenic_ranges$xmax, pathogenic_ranges$max_x_val, na.rm = T)

# Merge pathogenic ranges into main data
pathogenic_ranges = merge(pathogenic_ranges, data[,c("gene_sample", "gene")], all.x = T, by = 'gene')

#shape_cycle <- c(21, 22, 23, 24, 25, 16, 17, 15)

ggplot() +
  geom_tile(
    data = pathogenic_ranges,
    aes(
      x = (xmin + xmax) / 2,
      y = gene_sample,
      width  = xmax - xmin,
      height = 0.7,
      fill = "Pathogenic range"
    ),
    #fill = "grey85",
    color = NA
  ) + scale_fill_manual(
    name = NULL,
    values = c("Pathogenic range" = "grey85")
  ) +
  geom_jitter(
    data = data,
    aes(
      x = tool_allele_max,
      y = gene_sample,
      color = tool,
    ),
    size = 2,
    stroke = 1,
    width = 0, height = 0.3
  ) +
  geom_point(
    data = data,
    aes(
      x = molec_allele_max,
      y = gene_sample,
      shape = molec_type
    ),
    size = 2,
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Allele size (bp)",
    y = NULL,
    shape = "Molecular test",
    color = "Tool"
  ) +   
  facet_wrap(
    ~gene,
    scales = "free",
    space = "free_y",
    ncol = 1
  ) +
  scale_y_discrete(drop = TRUE) +
  scale_shape_discrete() +
  scale_shape_manual(
    values = rep(
      c(3, 4),
    ), na.translate = FALSE
  ) + 
  scale_color_manual(values = tool_colours)

ggsave('pathogenic_gene_strip.pdf', height = 12, width = 10) 



plot_data <- subset(data, gene %in% c('DMPK', 'FMR1', 'FXN', 'HTT', 'RFC1'))

facet_limits <- plot_data %>%
  group_by(gene) %>%
  # Find the absolute min and max considering both columns
  summarize(
    min_val = min(c(molec_allele_max, tool_allele_max), na.rm = TRUE),
    max_val = max(c(molec_allele_max, tool_allele_max), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Duplicate these values for both x and y mapping
  pivot_longer(cols = c(min_val, max_val), values_to = "range_val")

ggplot(plot_data, 
       aes(x = molec_allele_max, y = tool_allele_max, 
                 colour = tool, shape = molec_type)) + 
  geom_blank(data = facet_limits, aes(x = range_val, y = range_val), inherit.aes = FALSE) +
  
  geom_point() + geom_abline() + 
  scale_color_manual(values = tool_colours) +
  labs(x = "Molecular allele size (bp)", y = "Tool allele size (bp)",
       shape = "Molecular test", color = "Tool") +
  theme(aspect.ratio = 1) +
  facet_wrap(~gene, scales = "free")

ggsave('pathogenic_scatter.pdf', height = 7, width = 12) 





