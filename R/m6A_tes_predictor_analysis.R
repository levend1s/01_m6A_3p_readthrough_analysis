library(ggplot2)
library(scales)
library(dplyr)
library(ggExtra)
library(patchwork)
library(ggpubr)



args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
timepoint <- args[2]
out_dir <- args[3]

path_all <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/single_molecule_m6a_analysis.tsv"
timepoint <- "28"
out_dir <- "lala"

p_breaks <- c(0.01, 0.1, 0.5, 0.95, 0.99)

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

ggplot() +
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
    y = "Density"
  )

# ----- histogram: control vs knocksideways averages ----- #

ggplot() +
  geom_density(
    aes(x = group_ks_average$weighted_avg_probability_with_mods, fill = "methylated transcripts (knock-sideways)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_ks_average$weighted_avg_probability_without_mods, fill = "unmethylated transcripts (knock-sideways)"),
    alpha = alpha,
    linewidth = 1
  ) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_fill_manual(
  #   name   = "",
  #   values = c(
  #     "methylated transcripts (control)"   = "#1f77b4",
  #     "unmethylated transcripts (control)" = "#ff7f0e",
  #     "methylated transcripts (knock-sideways)" = "#1f77b4",
  #     "unmethylated transcripts (knock-sideways)" = "#ff7f0e"
  #   )
  # ) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  theme_classic(base_size = 20) +
  labs(
    x = "Probability",
    y = "Density"
  )


# ----- histogram: control vs knocksideways averages ----- #

ggplot() +
  geom_density(
    aes(x = group_control_average$weighted_avg_probability_with_mods, fill = "methylated transcripts (control)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_control_average$weighted_avg_probability_without_mods, fill = "unmethylated transcripts (control)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_ks_average$weighted_avg_probability_with_mods, fill = "methylated transcripts (knock-sideways)"),
    alpha = alpha,
    linewidth = 1
  ) +
  geom_density(
    aes(x = group_ks_average$weighted_avg_probability_without_mods, fill = "unmethylated transcripts (knock-sideways)"),
    alpha = alpha,
    linewidth = 1
  ) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_fill_manual(
  #   name   = "",
  #   values = c(
  #     "methylated transcripts (control)"   = "#1f77b4",
  #     "unmethylated transcripts (control)" = "#ff7f0e",
  #     "methylated transcripts (knock-sideways)" = "#1f77b4",
  #     "unmethylated transcripts (knock-sideways)" = "#ff7f0e"
  #   )
  # ) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  theme_classic(base_size = 20) +
  labs(
    x = "Probability",
    y = "Density"
  )


wilcox.test(c1_28_common$probability_with_mods, c2_28_common$probability_with_mods)
wilcox.test(k1_28_common$probability_with_mods, k2_28_common$probability_with_mods)
wilcox.test(c1_28_common$probability_with_mods, k1_28_common$probability_with_mods)

wilcox.test(group_control_average$weighted_avg_probability_with_mods, group_ks_average$weighted_avg_probability_with_mods)
wilcox.test(group_control_average$weighted_avg_probability_without_mods, group_ks_average$weighted_avg_probability_without_mods)

  