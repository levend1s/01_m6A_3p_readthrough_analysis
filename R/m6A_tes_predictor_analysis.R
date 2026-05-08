library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(ggsignif)

pval_to_signif <- function(p) {
  if (p < 0.0001) {
    return("****")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")   # not significant
  }
}

args <- commandArgs(trailingOnly = TRUE)
path_all <- args[1]
timepoint <- args[2]
ext = args[3]
out_dir <- args[4]

p_breaks <- c(0, 0.25, 0.5, 0.75, 1.0)

df_all <- read.delim(path_all, header=TRUE)
df_all$ID_start <- paste(df_all$ID, df_all$mod_start, sep = "-")
df_list <- split(df_all, df_all$sample)

# 28hpi
c1_28 <- df_list[[paste(timepoint, "C1_read_depth", sep="")]]
c2_28 <- df_list[[paste(timepoint, "C2_read_depth", sep="")]]
k1_28 <- df_list[[paste(timepoint, "K1_read_depth", sep="")]]
k2_28 <- df_list[[paste(timepoint, "K2_read_depth", sep="")]]

c1_28 <- c1_28[c1_28$percent_mod > 0.5,]
c2_28 <- c2_28[c2_28$percent_mod > 0.5,]

common_ids <- Reduce(intersect, list(
  c1_28$ID_start,
  c2_28$ID_start,
  k1_28$ID_start,
  k2_28$ID_start
))

c1_28_common <- c1_28[c1_28$ID_start %in% common_ids, ]
c2_28_common <- c2_28[c2_28$ID_start %in% common_ids, ]
k1_28_common <- k1_28[k1_28$ID_start %in% common_ids, ]
k2_28_common <- k2_28[k2_28$ID_start %in% common_ids, ]

control_28_average <- merge(c1_28_common, c2_28_common, by = "ID_start", suffixes = c("_1", "_2"))
ks_28_average <- merge(k1_28_common, k2_28_common, by = "ID_start", suffixes = c("_1", "_2"))

control_28_average$weighted_avg_probability_with_mods <-
  (control_28_average$probability_with_mods_1 * control_28_average$total_tx_with_mod_1 + 
     control_28_average$probability_with_mods_2 * control_28_average$total_tx_with_mod_2) /
  (control_28_average$total_tx_with_mod_1 + control_28_average$total_tx_with_mod_2)

control_28_average$weighted_avg_probability_without_mods <-
  (control_28_average$probability_without_mods_1 * control_28_average$total_tx_without_mod_1 + 
     control_28_average$probability_without_mods_2 * control_28_average$total_tx_without_mod_2) /
  (control_28_average$total_tx_without_mod_1 + control_28_average$total_tx_without_mod_2)

ks_28_average$weighted_avg_probability_with_mods <-
  (ks_28_average$probability_with_mods_1 * ks_28_average$total_tx_with_mod_1 + 
     ks_28_average$probability_with_mods_2 * ks_28_average$total_tx_with_mod_2) /
  (ks_28_average$total_tx_with_mod_1 + ks_28_average$total_tx_with_mod_2)

ks_28_average$weighted_avg_probability_without_mods <-
  (ks_28_average$probability_without_mods_1 * ks_28_average$total_tx_without_mod_1 + 
     ks_28_average$probability_without_mods_2 * ks_28_average$total_tx_without_mod_2) /
  (ks_28_average$total_tx_without_mod_1 + ks_28_average$total_tx_without_mod_2)

# assign test groups

group_c1 <- c1_28_common
group_c2 <- c2_28_common
group_k1 <- k1_28_common
group_k2 <- k2_28_common

group_control_average <- control_28_average
group_ks_average <- ks_28_average

# Plot 3 graphs: all samples for a timepoint
# (methylated x unmethylated) x group: averages
# knock sideways: averages unmethylated vs unmethylated

# ----- histogram: control vs knocksideways all samples ----- #

alpha <- 0.3
plot_name <- paste(out_dir, "/m6A_tes_predictor_all_samples", timepoint, "hpi.", ext, sep="")

out_file <- paste0(out_dir, "/", timepoint, "_all_samples_m6a_tes_predictor.", ext)

p <- ggplot() +
  geom_density(
    aes(x = group_c1$probability_with_mods, fill = "methylated transcripts (c1)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_k1$probability_with_mods, fill = "methylated transcripts (k1)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_c2$probability_with_mods, fill = "methylated transcripts (c2)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_k2$probability_with_mods, fill = "methylated transcripts (k2)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_c1$probability_without_mods, fill = "unmethylated transcripts (c1)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_k1$probability_without_mods, fill = "unmethylated transcripts (k1)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_c2$probability_without_mods, fill = "unmethylated transcripts (c2)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_k2$probability_without_mods, fill = "unmethylated transcripts (k2)"),
    alpha = alpha,
    linewidth = 1
  ) +
  scale_fill_brewer(palette = "Set1") +
  # scale_fill_manual(
  #   name   = "",
  #   values = c(
  #     "methylated transcripts"   = "#1f77b4",
  #     "unmethylated transcripts" = "#ff7f0e"
  #   )
  # ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  theme_classic(base_size = 20) +
  labs(
    x = "Probability",
    y = "Density",
    fill = ""
  )

p <- p +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.title = element_blank()
  )

ggsave(
  filename = out_file,   # output filename
  plot = p,                  # ggplot object
  device = ext,            # EPS format
  width = 10,                 # width in inches
  height = 4.8,                # height in inches
  units = "in",              # units (inches, cm, mm)
  bg = "transparent",
  dpi = 300                  # resolution, optional (ignored by EPS)
)

# ----- histogram: control vs knocksideways averages ----- #

out_file <- paste0(out_dir, "/", timepoint, "_knocksideways_m6a_tes_predictor.", ext)

num_bins <- 100

alpha <- 0.6

p <- ggplot() +
  geom_histogram(
    aes(x = group_ks_average$weighted_avg_probability_with_mods, fill = "methylated"),
    alpha = alpha,
    linewidth = 1,
    bins = num_bins
  ) +
  geom_histogram(
    aes(x = group_ks_average$weighted_avg_probability_without_mods, fill = "unmethylated"),
    alpha = alpha,
    linewidth = 1,
    bins = num_bins
  ) +
  scale_fill_brewer(palette = "Set2") +
  # scale_x_continuous(
  #   labels = percent_format(accuracy = 1)
  # ) +
  theme_classic(base_size = 20) +
  labs(
    x = "Probability",
    y = "Count",
    fill = ""
  )

p

p <- p +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.title = element_blank()
  )

ggsave(
  filename = out_file,   # output filename
  plot = p,                  # ggplot object
  device = ext,            # EPS format
  width = 7.4,                 # width in inches
  height = 4.8,                # height in inches
  units = "in",              # units (inches, cm, mm)
  bg = "transparent",
  dpi = 300                  # resolution, optional (ignored by EPS)
)


# ----- histogram: control vs knocksideways averages ----- #

out_file <- paste0(out_dir, "/", timepoint, "_sample_average_m6a_tes_predictor.", ext)

cm_ksm <- wilcox.test(group_control_average$weighted_avg_probability_with_mods, group_ks_average$weighted_avg_probability_with_mods)
cm_cu <- wilcox.test(group_control_average$weighted_avg_probability_with_mods, group_control_average$weighted_avg_probability_without_mods)
ksm_ksu <- wilcox.test(group_ks_average$weighted_avg_probability_with_mods, group_ks_average$weighted_avg_probability_without_mods)
cu_ksu <- wilcox.test(group_control_average$weighted_avg_probability_without_mods, group_ks_average$weighted_avg_probability_without_mods)

df <- bind_rows(
  tibble(
    probability = group_control_average$weighted_avg_probability_with_mods,
    condition = "control",
    methylation = "methylated"
  ),
  tibble(
    probability = group_control_average$weighted_avg_probability_without_mods,
    condition = "control",
    methylation = "unmethylated"
  ),
  tibble(
    probability = group_ks_average$weighted_avg_probability_with_mods,
    condition = "knock-sideways",
    methylation = "methylated"
  ),
  tibble(
    probability = group_ks_average$weighted_avg_probability_without_mods,
    condition = "knock-sideways",
    methylation = "unmethylated"
  )
)

p <- ggplot(
  df,
  aes(
    x = condition,
    y = probability,
    fill = methylation
  )
) +
  geom_violin(
    position = position_dodge(width = 0.8),
    alpha = alpha,
    linewidth = 1,
    trim = TRUE
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(
    # labels = percent_format(accuracy = 1),
    breaks = p_breaks
  ) +
  labs(
    x = "",
    # y = expression("Probability tx ends within 100nt of canonical m"^6 * "A"),
    y = "Probability",
    fill = ""
  ) +
  geom_signif(
    xmin = c(0.8, 1.8),
    xmax = c(1.2, 2.2),
    y_position = c(1.05, 1.05),
    annotations = c(pval_to_signif(cm_cu$p.value), pval_to_signif(ksm_ksu$p.value)),
    tip_length = 0.01
  ) +
  geom_signif(
    xmin = c(0.8),
    xmax = c(1.8),
    y_position = 1.15,
    annotations = pval_to_signif(cm_ksm$p.value),
    tip_length = 0.01
  ) +
  geom_signif(
    xmin = c(1.2),
    xmax = c(2.2),
    y_position = 1.25,
    annotations = pval_to_signif(cu_ksu$p.value),
    tip_length = 0.01
  ) +
  theme_classic(base_size = 20)

p

p <- p +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.title = element_blank()
  )

ggsave(
  filename = out_file,   # output filename
  plot = p,                  # ggplot object
  device = ext,            # EPS format
  bg = "transparent",
  width = 7.4,                 # width in inches
  height = 4.8,                # height in inches
  units = "in",              # units (inches, cm, mm)
  dpi = 300                  # resolution, optional (ignored by EPS)
)

# Do methylated transcripts across all samples come from the same distributuion? 
kruskal.test(list(group_c1$probability_with_mods, group_c2$probability_with_mods, group_k1$probability_with_mods, group_k2$probability_with_mods))

# Do UNMETHYLATED transcripts across all samples come from the same distributuion? 
kruskal.test(list(group_c1$probability_without_mods, group_c2$probability_without_mods, group_k1$probability_without_mods, group_k2$probability_without_mods))

  