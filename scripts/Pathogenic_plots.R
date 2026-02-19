library(ggplot2)
library(tidyverse)
library(cowplot)
library(purrr)
theme_set(theme_cowplot()) # Sets the default for subsequent plots

#colors <- c("#7FB3D5", "#F5A623", "#7DCEA0", "#F1948A", "#B39DDB", "#F7DC6F", "#B3911F")
#tableau10 
colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948","#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")
tools <- c("medaka", "atarva", "longtr", "strkit", "vamos", "straglr", "strdust")
tool_labels <- setNames(c("Medaka Tandem", "ATaRVa", "LongTR", "STRkit", "vamos", "Straglr", "STRdust"), tools)
tools <- sort(tools)
tool_colors <- setNames(colors, tools)

#data <- read_tsv("pathogenic_results.default_settings.tsv")
data <- read_tsv("pathogenic_results.tsv")

# Offset PABPN1 molecular size +4 to match STRchive coordinates
data[data$gene == 'PABPN1','molec_allele1_len'] = 4*3 + data[data$gene == 'PABPN1','molec_allele1_len']
data[data$gene == 'PABPN1','molec_allele2_len'] = 4*3 + data[data$gene == 'PABPN1','molec_allele2_len']

data$tool_allele_max = pmax(data$allele1_len, data$allele2_len, na.rm = T)
data$molec_allele_max = pmax(data$molec_allele1_len, data$molec_allele2_len, na.rm = T)
data$allele_max_diff = data$tool_allele_max - data$molec_allele_max
data$gene_sample = factor(interaction(data$gene, data$sample, sep = " | ", drop = TRUE))
data$molec_type[data$molec_type == "Unknown"] = NA

# data.long <- data %>%
#   rowwise() %>%
#   mutate(
#     molec_min = min(molec_allele1_len, molec_allele2_len, na.rm = TRUE),
#     molec_max = max(molec_allele1_len, molec_allele2_len, na.rm = TRUE),
#     allele_min = min(allele1_len, allele2_len, na.rm = TRUE),
#     allele_max = max(allele1_len, allele2_len, na.rm = TRUE)
#   ) %>%
#   ungroup()
#   
#   data %>%
#   rowwise() %>%
#   mutate(
#     molec_allele_len = list(sort(
#       c(molec_allele1_len, molec_allele2_len),
#       na.last = TRUE
#     )),
#     allele_len = list(sort(
#       c(allele1_len, allele2_len),
#       na.last = TRUE
#     ))
#   ) %>%
#   ungroup() %>% 
#   select(
#     -molec_allele1_len, -molec_allele2_len,
#     -allele1_len, -allele2_len
#   )

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
    values = c("Pathogenic range" = "grey90")
  ) +
  geom_point(
    data = data,
    aes(x = molec_allele1_len, y = gene_sample,
        shape = molec_type),
    size = 2,
  ) +
  geom_point(
    data = data,
    aes(x = molec_allele2_len, y = gene_sample,
        shape = molec_type),
    size = 2,
  ) +
  geom_jitter(
    data = data,
    aes(x = allele1_len, y = gene_sample, color = tool),
    size = 2, width = 0.3, height = 0.5, alpha = 0.8
  ) +
  geom_jitter(
    data = data,
    aes(x = allele2_len, y = gene_sample, color = tool),
    size = 2, width = 0.3, height = 0.5, alpha = 0.8
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
      c(15, 17),
    ), na.translate = FALSE
  ) + 
  scale_color_manual(
    values = tool_colors,
    labels = tool_labels
  )

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
  labs(x = "Molecular test, largest allele size (bp)", y = "Tool, largest allele size (bp)",
       shape = "Molecular test", color = "Tool") +
  theme(aspect.ratio = 1) +
  scale_shape_manual(
    values = rep(
      c(15, 17),
    ), na.translate = FALSE
  ) +
  scale_color_manual(
    values = tool_colors,
    labels = tool_labels
  ) +
  facet_wrap(~gene, scales = "free") 

ggsave('pathogenic_scatter.pdf', height = 7, width = 12) 


## Calculate sensitivity (here defined as identifying at least one pathogenic-sized allele in affected sample)
# Subset data to just those samples with pathogenic expansions
data_path = merge(data, unique(pathogenic_ranges[,c("gene", "path_min")]), all.x = T, by = 'gene')
data_path = subset(data_path, molec_allele_max >= path_min)[,c("gene", "tool", "sample", "path_min", 
                                                   "molec_allele1_len", "molec_allele2_len",
                                                   "allele1_len", "allele2_len", 
                                                   "tool_allele_max")]
affected_samples_count = length(unique(data_path$sample))
data_path$TP = data_path$tool_allele_max >= data_path$path_min
group_by(data_path, tool) %>% summarise(sum(TP)/affected_samples_count)

# Define set of preferred motifs
get_canonical <- function(x, preferred = c("CAG", "CTG", "CGG", "AAGGG", "GGGGCC")) {
  if (is.na(x) || x == "") return(x)
  
  # Generate all circular rotations
  n <- nchar(x)
  s_dup <- paste0(x, x)
  rotations <- sapply(1:n, function(i) substr(s_dup, i, i + n - 1))
  
  # 2. Check if any rotation is in your preferred list
  match <- intersect(rotations, preferred)
  
  if (length(match) > 0) {
    return(match[1]) # Return the first matching preferred motif found
  } else {
    return(min(rotations)) # Fallback to alphabetical min if no match
  }
}

# Count complete non-overlapping motifs of length n
count_motifs <- function(sequence, n, trim_incomplete = FALSE) {
  
  if (is.na(sequence) || is.na(n)) {
    return(NA)
  }
  
  sequence <- gsub("\\s+", "", sequence)
  seq_len <- nchar(sequence)
  
  if (n > seq_len) {
    return(NA)
  }
  
  # Trim trailing incomplete motif if requested
  if (trim_incomplete) {
    usable_len <- floor(seq_len / n) * n
    sequence <- substr(sequence, 1, usable_len)
    seq_len <- nchar(sequence)
  }
  
  # Extract non-overlapping motifs
  starts <- seq(1, seq_len - n + 1, by = n)
  motifs <- substring(sequence, starts, starts + n - 1)
  
  # Run-length encoding
  r <- rle(motifs)
  
  data.frame(
    motif = r$values,
    run_length = r$lengths,
    stringsAsFactors = FALSE
  )
}

count_motifs <- function(sequence, n, min_repeats = 2) {
  
  if (is.na(sequence) || is.na(n)) {
    return(NA)
  }
  
  sequence <- gsub("\\s+", "", sequence)
  seq_len <- nchar(sequence)
  
  if (n > seq_len) {
    return(NA)
  }
  
  # Step 1: count all sliding-window motifs
  starts <- 1:(seq_len - n + 1)
  motifs <- substring(sequence, starts, starts + n - 1)
  motif_counts <- sort(table(motifs), decreasing = TRUE)
  
  # Keep motifs that occur at least min_repeats times
  valid_motifs <- names(motif_counts[motif_counts >= min_repeats])
  
  if (length(valid_motifs) == 0) {
    return(data.frame(motif = sequence, run_length = 1))
  }
  
  results <- list()
  i <- 1
  
  while (i <= seq_len) {
    
    matched <- FALSE
    
    # Try matching any valid motif (in order of frequency)
    if (i <= seq_len - n + 1) {
      for (motif in valid_motifs) {
        if (substr(sequence, i, i + n - 1) == motif) {
          
          run_len <- 0
          
          while (i <= seq_len - n + 1 &&
                 substr(sequence, i, i + n - 1) == motif) {
            run_len <- run_len + 1
            i <- i + n
          }
          
          results[[length(results) + 1]] <-
            data.frame(motif = motif,
                       run_length = run_len,
                       stringsAsFactors = FALSE)
          
          matched <- TRUE
          break
        }
      }
    }
    
    # If no valid motif matched → accumulate interruption
    if (!matched) {
      start_int <- i
      
      while (i <= seq_len &&
             !(i <= seq_len - n + 1 &&
               any(substr(sequence, i, i + n - 1) == valid_motifs))) {
        i <- i + 1
      }
      
      interruption_seq <- substr(sequence, start_int, i - 1)
      
      results[[length(results) + 1]] <-
        data.frame(motif = interruption_seq,
                   run_length = 1,
                   stringsAsFactors = FALSE)
    }
  }
  
  t = do.call(rbind, results)
  #return(summarize_repeat_structure(do.call(rbind, results)))
  
  #return(arrange(t, desc(run_length))[1,])
  l = list(names(motif_counts[1]), motif_counts[1])
  
  #give names to elements of list
  names(l) <- c('motif', 'count') 
  return(l)
}

summarize_repeat_structure <- function(run_df) {
  
  if (is.null(run_df) || all(is.na(run_df))) {
    return(NA)
  }
  
  if (!all(c("motif", "run_length") %in% colnames(run_df))) {
    stop("Input must have columns 'motif' and 'run_length'")
  }
  
  # Determine repeat motif length (assume most common length > 1)
  motif_lengths <- nchar(run_df$motif)
  n <- as.numeric(names(sort(table(motif_lengths[motif_lengths > 1]),
                             decreasing = TRUE))[1])
  
  # Separate repeat motifs from interruptions
  is_repeat <- nchar(run_df$motif) == n
  
  repeat_df <- run_df[is_repeat, ]
  interruption_df <- run_df[!is_repeat, ]
  
  # Sum total occurrences per motif
  repeat_counts <- aggregate(run_length ~ motif,
                             data = repeat_df,
                             sum)
  
  # Count interruption blocks
  n_interruptions <- nrow(interruption_df)
  
  return(n_interruptions)
  
  # list(
  #   motif_counts = repeat_counts,
  #   n_interruptions = n_interruptions
  # )
}



data = merge(data, unique(pathogenic_ranges[,c("gene", "motif_size")]), all.x = T, by = 'gene')
data <- data %>%
  mutate(
    allele1_main = map2(allele1_seq, motif_size,
                                    ~ count_motifs(.x, .y)),
    allele2_main = map2(allele2_seq, motif_size,
                                    ~ count_motifs(.x, .y))
  ) %>%
      # Splits the named list into columns with a prefix
      unnest_wider(allele1_main, names_sep = "_") %>%
      unnest_wider(allele2_main, names_sep = "_") %>%
  select(-ends_with("main_1"))

data <- data %>%
  mutate(
    allele1_main_prop = (allele1_main_count * motif_size)/allele1_len,
    allele2_main_prop = (allele2_main_count * motif_size)/allele2_len
  )

data <- data %>%
  mutate(
    allele1_main_motif = map_chr(allele1_main_motif, get_canonical),
    allele1_main_motif = map_chr(allele2_main_motif, get_canonical)
  )





  # ) %>%
  # mutate(allele1_interruptions = map(allele1_motifs_counts,
  #                               ~ summarize_repeat_structure(.x)),
  #        allele2_interruptions = map(allele2_motifs_counts,
  #                               ~ summarize_repeat_structure(.x))
  #        
  #        )

ggplot(data) + 
  geom_jitter(aes(y = gene_sample, x = allele1_main_prop,
                 color = tool),
              width = 0.01, height = 0.3, alpha = 0.8) +
  geom_jitter(aes(y = gene_sample, x = allele2_main_prop,
                 color = tool),
              width = 0.01, height = 0.3, alpha = 0.8) +
  scale_color_manual(
    values = tool_colors,
    labels = tool_labels
  ) +   
  facet_wrap(
    ~gene,
    scales = "free",
    space = "free_y",
    ncol = 1
  ) +
  scale_y_discrete(drop = TRUE) 

ggsave('pathogenic_main_motif.pdf', height = 14, width = 10) 

ggplot(data) + 
  geom_jitter(aes(color = allele1_main_motif, x = allele1_main_prop, y = gene_sample)) + 
  geom_jitter(aes(color = allele2_main_motif, x = allele2_main_prop, y = gene_sample)) + 
  facet_wrap(
  ~gene,
  scales = "free",
  ) +
  scale_y_discrete(drop = TRUE, na.translate = FALSE) 
