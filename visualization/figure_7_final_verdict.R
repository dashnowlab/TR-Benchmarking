# ============================================================
# HEATMAP ONLY (no winners, no extra comments)
# Fixes:
#  - Category order: Usability first, then Accuracy
#  - Metric order in Accuracy: Assembly (3 rows) -> Mendelian -> Pathogenic
#  - Legend shows 1.00 at top
#  - Save to both PNG and PDF
#  - Heatmap text increased by 1 point
#  - Legend spacing fixed
# ============================================================

library(tidyverse)
library(scales)
library(grid)

df <- tribble(
  ~Category,   ~Metric,                                  ~ATaRVa, ~LongTR, ~`Medaka Tandem`, ~Straglr, ~STRdust, ~STRkit, ~vamos,
  "Usability", "Mean runtime (mins)",                       199,      98,        209,          813,      234,      90,      73,
  "Usability", "Max memory use (GB)",                       1.5,      6.0,       3.1,          70.0,     1.6,     3.1,     1.3,
  "Accuracy",  "Assembly Consistency \n ≤ 1 motif diff (%)",   99.41,    99.38,     99.37,        93.66,    94.93,   98.47,   99.17,
  "Accuracy",  "Assembly Consistency \n Homopolymers (%)",     79.38,    79.19,     89.36,        11.39,    31.00,   58.19,   76.49,
  "Accuracy",  "Assembly Sensitivity \n Expansions >500 bp (%)",85,       87,        75,           72,       77,      35,      83,
  "Accuracy",  "Mendelian Consistency \n ≤ 1 motif diff (%)",  99.01,    98.01,     99.17,        99.20,    74.46,   97.71,   97.23
)

df_path <- tribble(
  ~Category,  ~Metric,                           ~ATaRVa, ~LongTR, ~`Medaka Tandem`, ~Straglr, ~STRdust, ~STRkit, ~vamos,
  "Accuracy", "Pathogenic loci \n sensitivity (%)", 70,      87,        39,           61,       96,      61,      61
)

df <- bind_rows(df, df_path)

tools_order <- c("ATaRVa","LongTR","Medaka Tandem","Straglr","STRdust","STRkit","vamos")
lower_is_better <- c("Mean runtime (mins)", "Max memory use (GB)")

tableau10 <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
)

t_red   <- tableau10[6]
t_green <- tableau10[5]

stop_vals <- c(0.00, 0.25, 0.50, 1.00)
stop_cols <- c(t_red, t_red, "white", t_green)

legend_breaks <- c(1e-3, 0.25, 0.50, 0.75, 1 - 1e-3)
legend_labels <- c("0.00", "0.25", "0.50", "0.75", "1.00")

metric_levels <- c(
  "Max memory use (GB)",
  "Mean runtime (mins)",
  "Pathogenic loci \n sensitivity (%)",
  "Mendelian Consistency \n ≤ 1 motif diff (%)",
  "Assembly Sensitivity \n Expansions >500 bp (%)",
  "Assembly Consistency \n Homopolymers (%)",
  "Assembly Consistency \n ≤ 1 motif diff (%)"
)

dat <- df %>%
  pivot_longer(cols = all_of(tools_order), names_to="Tool", values_to="Value") %>%
  group_by(Metric) %>%
  mutate(
    norm  = rescale(Value, to=c(0,1), from=range(Value, na.rm=TRUE)),
    Score = if_else(Metric %in% lower_is_better, 1 - norm, norm),
    Label = case_when(
      str_detect(Metric, "\\(mins\\)") ~ sprintf("%.0f", Value),
      str_detect(Metric, "\\(GB\\)")   ~ sprintf("%.1f", Value),
      TRUE                             ~ sprintf("%.2f", Value)
    )
  ) %>%
  ungroup() %>%
  mutate(
    Tool     = factor(Tool, levels = tools_order),
    Category = factor(Category, levels = c("Usability","Accuracy")),
    Metric   = factor(Metric, levels = metric_levels),
    Score_fill = pmin(pmax(Score, 1e-3), 1 - 1e-3),
    Label = if_else(Metric == "Pathogenic loci \n sensitivity (%)" & Tool == "vamos",
                    "30–61*", Label)
  )

p <- ggplot(dat, aes(Tool, Metric, fill = Score_fill)) +
  geom_tile(color="white", linewidth=0.1) +
  geom_text(aes(label = Label), size=4.0, fontface="bold") +
  facet_grid(Category ~ ., scales="free_y", space="free_y") +
  scale_fill_gradientn(
    colours = stop_cols,
    values  = scales::rescale(stop_vals),
    trans   = scales::transform_logit(),
    limits  = c(1e-3, 1 - 1e-3),
    breaks  = legend_breaks,
    labels  = legend_labels,
    oob     = scales::squish,
    name    = "Relative\nperformance",
    guide   = guide_colorbar(
      direction = "vertical",
      barheight = unit(6, "cm"),
      barwidth  = unit(0.8, "cm"),
      ticks.colour = "black",
      frame.colour = NA
    )
  ) +
  scale_x_discrete(expand = expansion(mult=c(0.01,0.01))) +
  scale_y_discrete(expand = expansion(mult=c(0.01,0.01))) +
  labs(x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face="bold", size=10),
    axis.text.y = element_text(face="bold", size=10),
    strip.text.y = element_text(face="bold", size=10),
    strip.background = element_rect(fill="grey95", color=NA),
    panel.spacing.y = unit(0.08, "lines"),
    plot.margin = margin(4, 6, 4, 4),
    legend.title = element_text(face="bold", size=10),
    legend.text  = element_text(face="bold", size=9),
    legend.key.height = unit(1.2, "cm")
  )

print(p)

ggsave("tool_summary_heatmap_tableau_high_contrast.png", p, width=9.5, height=4.5, dpi=300)
ggsave("tool_summary_heatmap_tableau_high_contrast.pdf", p, width=9.5, height=4.5, device="pdf")