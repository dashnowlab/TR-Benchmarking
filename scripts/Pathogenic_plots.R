library(ggplot2)
library(tidyverse)
library(cowplot)
library(purrr)
theme_set(theme_cowplot()) # Sets the default for subsequent plots

# Set consistent tool colors
tableau_colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#FF9DA7")
tools <- c("longtr", "atarva", "medaka", "strkit", "vamos", "strdust", "straglr")
labels <- c("LongTR", "ATaRVa", "Medaka\nTandem", "STRkit", "vamos", "STRdust", "Straglr")

tool_levels <- tools
tool_labels <- setNames(labels, tools)
tool_cols <- setNames(tableau_colors, tools) # <-- map colors by tool KEYS

# analysis = 'minreads'
data <- read_tsv("pathogenic_results.min_reads.tsv")

# # analysis = 'default'
default.data <- read_tsv("pathogenic_results.default.tsv")
# straglr.default = subset(default.data, tool == 'straglr')
# 
# # Replace minreads straglr with default straglr
# data = subset(data, tool != 'straglr')
# data = rbind(data, straglr.default)

# analysis = 'best'
# old.data <- read_tsv("old/pathogenic_results.backup.tsv")
# vamos.old = subset(old.data, tool == 'vamos')

# # Replace new vamos with old vamos
# data = subset(data, tool != 'vamos')
# data = rbind(data, vamos.old)

prep_data <- function(mydata) {
  # Offset PABPN1 molecular size +4 to match STRchive coordinates
  mydata[mydata$gene == 'PABPN1','molec_allele1_len'] = 4*3 + mydata[mydata$gene == 'PABPN1','molec_allele1_len']
  mydata[mydata$gene == 'PABPN1','molec_allele2_len'] = 4*3 + mydata[mydata$gene == 'PABPN1','molec_allele2_len']
  
  mydata$tool_allele_max = pmax(mydata$allele1_len, mydata$allele2_len, na.rm = T)
  mydata$tool_allele_min = pmin(mydata$allele1_len, mydata$allele2_len, na.rm = T)
  mydata$molec_allele_max = pmax(mydata$molec_allele1_len, mydata$molec_allele2_len, na.rm = T)
  mydata$allele_max_diff = mydata$tool_allele_max - mydata$molec_allele_max
  mydata$gene_sample = factor(interaction(mydata$gene, mydata$sample, sep = " | ", drop = TRUE))
  mydata$molec_type[mydata$molec_type == "Unknown"] = NA
  
  return(mydata)
}

data = prep_data(data)
default.data = prep_data(default.data)

################################
# Pathogenic loci "strip" plot
################################

# Positions of shaded rectangles
pathogenic_ranges <- data.frame(
  gene = c(
          "HTT", "RFC1", "ATXN1", 
          "C9orf72", "DMPK", 
          "FMR1", "FMR1",
          "FXN", "NOTCH2NLC", "PABPN1"),
  path_min = c(
           36, 400, 39,
           31, 50, 
           45, 201,
           56, 66, 12),
  path_max = c(
           250, 2750, 91,
           4088, 4000, 
           200, 2000,
           1700, 517, 18),
  motif_size = c(3, 5, 3,
                 6, 3, 
                 3, 3,
                 3, 3, 3),
  size_range = c("Pathogenic", "Pathogenic", "Pathogenic",
                 "Pathogenic", "Pathogenic", 
                 "Premutation", "Pathogenic",
                 "Pathogenic", "Pathogenic", "Pathogenic")
)
pathogenic_ranges$xmin = pathogenic_ranges$path_min #* pathogenic_ranges$motif_size
pathogenic_ranges$xmax = pathogenic_ranges$path_max #* pathogenic_ranges$motif_size

# Add alleles in motif counts
strip_data = merge(data, unique(pathogenic_ranges[c("gene", "motif_size")]), all.x = T, by = 'gene')

# Divide all alleles by motif size to get into the same context as pathogenic thresholds
allele_cols = c("molec_allele1_len","molec_allele2_len","allele1_len","allele2_len","tool_allele_max","tool_allele_min", "molec_allele_max")
strip_data[,allele_cols] = round(strip_data[,allele_cols]/strip_data$motif_size)

gene_max_x <- strip_data %>%
  summarise(
    max_x_val = max(tool_allele_max, molec_allele_max, na.rm = TRUE),
    .by = gene
  )

# Set max point on plot
# max_x_val = max(data$tool_allele_max, data$molec_allele_max, na.rm = T)
pathogenic_ranges = merge(pathogenic_ranges, gene_max_x)
pathogenic_ranges$xmax = pmin(pathogenic_ranges$xmax, pathogenic_ranges$max_x_val, na.rm = T)

# Merge pathogenic ranges into main data
pathogenic_ranges = merge(pathogenic_ranges, strip_data[,c("gene_sample", "gene")], all.x = T, by = 'gene')

strip_data = unique(strip_data)
pathogenic_ranges = unique(pathogenic_ranges)

strip.plot = ggplot() +
  geom_tile(
    data = pathogenic_ranges,
    aes(
      x = (xmin + xmax) / 2,
      y = gene_sample,
      width  = xmax - xmin,
      height = 0.7,
      fill = size_range
    ),
    #fill = "grey85",
    color = NA
  ) + scale_fill_manual(
    name = "Allele size range",
    values = c("Premutation" = "grey95", "Pathogenic" = "grey85")
  ) +
  geom_point(
    data = strip_data,
    aes(x = molec_allele1_len, y = gene_sample,
        shape = molec_type),
    size = 2,
  ) +
  geom_point(
    data = strip_data,
    aes(x = molec_allele2_len, y = gene_sample,
        shape = molec_type),
    size = 2,
  ) +
  geom_jitter(
    data = strip_data,
    aes(x = allele1_len, y = gene_sample, color = tool),
    size = 2, alpha = 0.8#, width = 0.3, height = 0.5
  ) +
  geom_jitter(
    data = strip_data,
    aes(x = allele2_len, y = gene_sample, color = tool),
    size = 2, alpha = 0.8#, width = 0.3, height = 0.5
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position=c(-0.21,0.8),
    strip.text = element_blank()
  ) +
  labs(
    x = "Allele size (motifs)",
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
      c(15, 17),
    ), na.translate = FALSE
  ) + 
  scale_color_manual(
    values = tool_cols,
    breaks = tool_levels,
    labels = tool_labels[tool_levels],
    name = "Tool",
    drop = FALSE
  )

#ggsave(paste0('pathogenic_gene_strip.', analysis, '.pdf'), height = 12, width = 10) 

################################
# Calculate sensitivity
################################
# Sensitivity is here defined as identifying at least one pathogenic-sized allele
# in heterozygous samples and two in homozygous expanded samples

# # Build new df to subset to just those samples with pathogenic expansions
# sensitivity_path_ranges <- unique(pathogenic_ranges[,c("gene", "path_min", "motif_size")]) %>%
#   group_by(gene) %>%
#   slice_max(path_min, with_ties = FALSE) %>%
#   ungroup()
# 
# data_path = merge(data, sensitivity_path_ranges, all.x = T, by = 'gene')
# # Divide all alleles by motif size to get into the same context as pathogenic thresholds
# allele_cols = c("molec_allele1_len","molec_allele2_len","allele1_len","allele2_len","tool_allele_max","tool_allele_min")
# data_path[,allele_cols] = round(data_path[,allele_cols]/data_path$motif_size)
# # Different pathogenic threshold for premutation samples
# premutation_samples = c("NA06905", "CD00014")
# data_path$path_min[data_path$sample %in% premutation_samples] = unique(pathogenic_ranges$path_min[pathogenic_ranges$size_range == "Premutation" & pathogenic_ranges$gene == "FMR1"])
# 
# # Remove samples with benign alleles
# #data_path = subset(data_path, molec_allele_max >= path_min)
# controls = c("R210029", "NA06890", "NA06893", "CD00022")
# data_path = data_path[!(data_path$sample %in% controls),]
# 
# # Simplify data frame
# data_path = data_path[,c("gene", "tool", "sample",  "gene_sample", 
#                                                                "motif_size", "path_min", 
#                                                                "molec_allele1_len", "molec_allele2_len",
#                                                                "allele1_len", "allele2_len", 
#                                                                "tool_allele_max", "tool_allele_min")]
# 
# 
# # Count samples with expansions or premutations
# affected_samples_count = length(unique(data_path$sample))
# 
# data_path$TP = data_path$tool_allele_max >= data_path$path_min
# # homozygous samples
# hom_samples = c("NA15850","R190082","R210005","R210020","R210023")
# data_path$TP[data_path$sample %in% hom_samples] = (data_path$tool_allele_min >= data_path$path_min)[data_path$sample %in% hom_samples]
# 
# # Calculate per-sample TPs
# group_by(data_path, tool, sample) %>% slice_max(tool, with_ties = FALSE) %>% ungroup() -> per_sample_TPs
# 
# ## Write overall sensitivity to table
# group_by(per_sample_TPs,tool) %>%
#   summarise(sensitivity = round(sum(TP)/affected_samples_count, digits = 2)) -> sensitivity

calc_TPs <- function(in.data, pathogenic_ranges) {
  # Build new df to subset to just those samples with pathogenic expansions
  sensitivity_path_ranges <- unique(pathogenic_ranges[,c("gene", "path_min", "motif_size")]) %>%
    group_by(gene) %>%
    slice_max(path_min, with_ties = FALSE) %>%
    ungroup()
  
  data_path = merge(in.data, sensitivity_path_ranges, all.x = T, by = 'gene')
  # Divide all alleles by motif size to get into the same context as pathogenic thresholds
  allele_cols = c("molec_allele1_len","molec_allele2_len","allele1_len","allele2_len","tool_allele_max","tool_allele_min")
  data_path[,allele_cols] = round(data_path[,allele_cols]/data_path$motif_size)
  # Different pathogenic threshold for premutation samples
  premutation_samples = c("NA06905", "CD00014")
  data_path$path_min[data_path$sample %in% premutation_samples] = unique(pathogenic_ranges$path_min[pathogenic_ranges$size_range == "Premutation" & pathogenic_ranges$gene == "FMR1"])
  
  # Remove samples with benign alleles
  #data_path = subset(data_path, molec_allele_max >= path_min)
  controls = c("R210029", "NA06890", "NA06893", "CD00022")
  data_path = data_path[!(data_path$sample %in% controls),]
  
  # Simplify data frame
  data_path = data_path[,c("gene", "tool", "sample",  "gene_sample", 
                           "motif_size", "path_min", 
                           "molec_allele1_len", "molec_allele2_len",
                           "allele1_len", "allele2_len", 
                           "tool_allele_max", "tool_allele_min")]
  
  
  # Count samples with expansions or premutations
  affected_samples_count = length(unique(data_path$sample))
  
  data_path$TP = data_path$tool_allele_max >= data_path$path_min
  # homozygous samples
  hom_samples = c("NA15850","R190082","R210005","R210020","R210023")
  data_path$TP[data_path$sample %in% hom_samples] = (data_path$tool_allele_min >= data_path$path_min)[data_path$sample %in% hom_samples]
  
  # Calculate per-sample TPs
  group_by(data_path, tool, sample) %>% slice_max(tool, with_ties = FALSE) %>% ungroup() -> per_sample_TPs
  
  return(per_sample_TPs)
}

calc_sensitivity <- function(in.data, pathogenic_ranges) {
  # Build new df to subset to just those samples with pathogenic expansions
  sensitivity_path_ranges <- unique(pathogenic_ranges[,c("gene", "path_min", "motif_size")]) %>%
    group_by(gene) %>%
    slice_max(path_min, with_ties = FALSE) %>%
    ungroup()
  
  data_path = merge(in.data, sensitivity_path_ranges, all.x = T, by = 'gene')
  # Divide all alleles by motif size to get into the same context as pathogenic thresholds
  allele_cols = c("molec_allele1_len","molec_allele2_len","allele1_len","allele2_len","tool_allele_max","tool_allele_min")
  data_path[,allele_cols] = round(data_path[,allele_cols]/data_path$motif_size)
  # Different pathogenic threshold for premutation samples
  premutation_samples = c("NA06905", "CD00014")
  data_path$path_min[data_path$sample %in% premutation_samples] = unique(pathogenic_ranges$path_min[pathogenic_ranges$size_range == "Premutation" & pathogenic_ranges$gene == "FMR1"])
  
  # Remove samples with benign alleles
  #data_path = subset(data_path, molec_allele_max >= path_min)
  controls = c("R210029", "NA06890", "NA06893", "CD00022")
  data_path = data_path[!(data_path$sample %in% controls),]
  
  # Simplify data frame
  data_path = data_path[,c("gene", "tool", "sample",  "gene_sample", 
                           "motif_size", "path_min", 
                           "molec_allele1_len", "molec_allele2_len",
                           "allele1_len", "allele2_len", 
                           "tool_allele_max", "tool_allele_min")]
  
  
  # Count samples with expansions or premutations
  affected_samples_count = length(unique(data_path$sample))
  
  data_path$TP = data_path$tool_allele_max >= data_path$path_min
  # homozygous samples
  hom_samples = c("NA15850","R190082","R210005","R210020","R210023")
  data_path$TP[data_path$sample %in% hom_samples] = (data_path$tool_allele_min >= data_path$path_min)[data_path$sample %in% hom_samples]
  
  # Calculate per-sample TPs
  group_by(data_path, tool, sample) %>% slice_max(tool, with_ties = FALSE) %>% ungroup() -> per_sample_TPs
  
    
  ## Write overall sensitivity to table
  group_by(per_sample_TPs,tool) %>%
  summarise(sensitivity = round(sum(TP)/affected_samples_count, digits = 2)) -> sensitivity
    
  return(sensitivity)
}


sensitivity.minreads = calc_sensitivity(data, pathogenic_ranges)
colnames(sensitivity.minreads) = c('tool', 'minreads')
sensitivity.minreads$tool[sensitivity.minreads$tool == 'vamos'] = 'vamos v3.0.5'
sensitivity.minreads$tool[sensitivity.minreads$tool == 'vamos_old'] = 'vamos v2.1.7'

sensitivity.default = calc_sensitivity(default.data, pathogenic_ranges)
colnames(sensitivity.default) = c('tool', 'default')
sensitivity.default$tool[sensitivity.default$tool == 'vamos'] = 'vamos v3.0.5'
sensitivity.default$tool[sensitivity.default$tool == 'vamos_old'] = 'vamos v2.1.7'

sensitivity.table = merge(sensitivity.default, sensitivity.minreads, all = T)
sensitivity.table$default[sensitivity.table$tool == 'vamos v2.1.7'] = sensitivity.table$minreads[sensitivity.table$tool == 'vamos v2.1.7']
sensitivity.table$minreads[sensitivity.table$tool %in% c('vamos v2.1.7', 'vamos v3.0.5')] = NA
sensitivity.table$default[sensitivity.table$tool == 'longtr']  = 0

sensitivity.table$minreads = scales::percent(sensitivity.table$minreads)
sensitivity.table$default = scales::percent(sensitivity.table$default)

write_tsv(sensitivity.table, 'Sensitivity.tsv')

# Plot the table
library(ggpubr)

sensitivity.table.plot = ggtexttable(sensitivity.table, rows = NULL, 
                                theme = ttheme("light", base_size = 16, padding = unit(c(5, 5), "mm")))


## Plot all the genotypes coloured by TP true/false
per_sample_TPs = calc_TPs(data, pathogenic_ranges)
# Set plotting order of everything
per_sample_TPs$gene_sample = factor(per_sample_TPs$gene_sample, levels = sort(levels(per_sample_TPs$gene_sample), decreasing = T)) 
# Create a molecular genotype for plotting
per_sample_TPs = unite(per_sample_TPs, "molec", molec_allele1_len, molec_allele2_len, sep = "/", na.rm = TRUE, remove = FALSE)
per_sample_TPs = unite(per_sample_TPs, "text_genotype", allele1_len, allele2_len, sep = "/", na.rm = TRUE, remove = FALSE)

plot_TPs = subset(per_sample_TPs, tool != "vamos")
plot_TPs$tool[plot_TPs$tool == 'vamos_old'] = 'vamos'
plot.sensitivity.minreads = sensitivity.minreads
plot.sensitivity.minreads$tool[plot.sensitivity.minreads$tool == 'vamos v2.1.7'] = 'vamos'
plot.sensitivity.minreads$minreads = scales::percent(plot.sensitivity.minreads$minreads)
x_axis_order = c("pathogenic\nmin", "molec", sort(unique(c(plot_TPs$tool, "longtr"))))

# Sets default geom_text size to approximately 12pt
update_geom_defaults("text", list(size = 12 / .pt))

sensitivity.plot = ggplot(plot_TPs, aes(x = tool, y = gene_sample)) +
  geom_tile(aes(fill = TP)) + 
  geom_tile(aes(y = 'Total sensitivity', fill = 'Total sensitivity')) +
  geom_text(aes(label = text_genotype)) + 
  geom_text(aes(x = "pathogenic\nmin", label = path_min)) +
  geom_text(aes(x = "molec", label = molec)) +
  geom_text(data = plot.sensitivity.minreads, aes(y = 'Total sensitivity', x = tool, label = minreads)) +
  labs(y = 'Gene | Sample', x = "Method (allele sizes in motifs)") + 
  scale_x_discrete(limits = x_axis_order) +
  scale_fill_manual(
    values = c("#98df8a", "#ff9896", "#aec7e8"),
    labels = c("True positive", "False negative", "Total sensitivity"),
    breaks = c("TRUE", "FALSE", "Total sensitivity")
  ) +
  theme(legend.position="top", legend.title = element_blank(),
        text = element_text(size = 14))

ggsave('Sensitivity.pdf', sensitivity.plot, height = 8, width = 12)

left.panel = plot_grid(sensitivity.table.plot, sensitivity.plot, 
                       labels = c("A", "B"), 
                       ncol = 1, nrow = 2, 
                       rel_heights = c(1, 1.7),
                       align = 'h', axis = 'l')
multipanel.fig <- plot_grid(left.panel, strip.plot, 
                            labels = c("", "C"), 
                            ncol = 2, nrow = 1, 
                            rel_widths = c(1, 1.2))
ggsave('Pathogenic_figure.pdf', plot = multipanel.fig, height = 12, width = 20)

