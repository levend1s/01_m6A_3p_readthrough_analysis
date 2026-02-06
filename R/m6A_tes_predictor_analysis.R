library(ggplot2)
library(scales)
library(dplyr)
library(ggExtra)
library(patchwork)
library(ggpubr)

path_ctrl <- "/Users/jlevendis/rqc/100nt_28C2_m6A_tes_predictor.tsv"
path_ks <- "/Users/jlevendis/rqc/100nt_28K1_m6A_tes_predictor.tsv"
# path <- "/Users/jlevendis/rqc/100nt0.65p_28K1_m6A_tes_predictor.tsv"
# path <- "/Users/jlevendis/rqc/m6A_tes_predictor.tsv"

# path <- "/Users/jlevendis/rqc/1p_100nt_28C1_m6A_tes_predictor.tsv"
# path <- "/Users/jlevendis/rqc/test.tsv"


df_ctrl <- read.delim(path_ctrl, header=TRUE)
df_ks <- read.delim(path_ks, header=TRUE)

# df <- df[df$percent_mod > 0.1,]
df_ctrl <- df_ctrl[df_ctrl$percent_mod > 0.5,]
df_ctrl$ID_start <- paste(df_ctrl$ID, df_ctrl$mod_start, sep = "-")
df_ks$ID_start <- paste(df_ks$ID, df_ks$mod_start, sep = "-")

df_common <- merge(df_ctrl, df_ks, by = c("ID_start"))
# df_ks <- df_ks[df_ks$percent_mod > 0.5,]
# df <- df[df$probability > 0.01 & df$probability < 0.99,]

p_breaks <- c(0.01, 0.1, 0.5, 0.95, 0.99)

# ----- histogram: percent modified ----- # 
# ggplot(df, aes(x = percent_mod)) +
#   geom_histogram(bins = 30, fill = "steelblue", color = "white") +
#   theme_classic(base_size = 20) +
#   scale_x_continuous(
#     labels  = percent_format(accuracy = 1)
#   ) +
#   labs(
#     x = "Percent modified",
#     y = "Count"
#   )



# ----- histogram: probability transcript containing an m6A ends within 100nt ----- #
# ggplot(df, aes(x = probability_with_mods)) +
#   geom_histogram(bins = 30, fill = "steelblue", color = "white") +
#   theme_classic(base_size = 20) +
#   scale_x_continuous(
#     labels  = percent_format(accuracy = 1)
#   ) +
#   labs(
#     x = "probability",
#     y = "Count"
#   )
# 
# num_bins = 30
# 
# ggplot(df_common) +
#   geom_histogram(
#     aes(x = probability_with_mods.x, fill = "methylated transcripts"),
#     bins = num_bins,
#     alpha = 0.5,
#     position = "identity"
#   ) +
#   geom_histogram(
#     aes(x = probability_without_mods.x, fill = "unmethylated transcripts"),
#     bins = num_bins,
#     alpha = 0.5,
#     position = "identity"
#   ) +
#   scale_fill_manual(
#     name   = "",
#     values = c("methylated transcripts" = "#1f77b4", "unmethylated transcripts" = "#ff7f0e")
#   ) +
#   scale_x_continuous(
#     labels  = percent_format(accuracy = 1)
#   ) +
#   # scale_y_continuous(
#   #   transform = "log2"
#   # ) +
#   theme_classic(base_size = 20) +
#   labs(
#     x = "Probability",
#     y = "Count",
#     title = "Do transcripts end within 100nt of m6A sites?"
#   )


# ----- scatter: num_tx_within vs probability ----- #
# ggplot() +
#   geom_point(data = df, aes(x = probability_with_mods, y = probability_without_mods), alpha = 1, size=1) +
#   scale_y_continuous(
#     # trans = logit_trans(),
#     # breaks  = p_breaks,
#     labels  = percent_format(accuracy = 1),
#     name = "Probability (unmethylated transcripts)"
#   ) +
#   scale_x_continuous(
#     # trans = logit_trans(),
#     # breaks  = p_breaks,
#     labels  = percent_format(accuracy = 1),
#     name = "Probability (methylated transcripts)"
#   ) +
#   theme_classic(base_size = 20)

# # ----- 2d heatmap: probability vs percent mod (raw) ----- #
# ggplot(df, aes(x = probability_without_mods, y = percent_mod)) +
#   geom_bin2d(bins = 100) +
#   scale_x_continuous(
#     labels  = percent_format(accuracy = 1)
#   ) +
#   scale_y_continuous(
#     labels  = percent_format(accuracy = 1)
#   ) +
#   scale_fill_viridis_c(
#     trans = "log10",
#     name = "num sites"
#   ) +
#   labs(
#     x = "Probability transcripts end within 200nt of m6A",
#     y = "m6A / A"
#   ) +
#   theme_classic(base_size = 20)

# ----- 2d heatmap: probability with vs without mods (raw) ----- #
# ggplot(df, aes(x = probability_without_mods, y = percent_mod)) +
#   geom_bin2d(bins = 30) +
#   scale_x_continuous(
#     labels  = percent_format(accuracy = 1)
#   ) +
#   scale_y_continuous(
#     labels  = percent_format(accuracy = 1)
#   ) +
#   scale_fill_viridis_c(
#     trans = "log10",
#     name = "num sites"
#   ) +
#   labs(
#     x = "Probability (unmethylated)",
#     y = "Probability (methylated)"
#   ) +
#   theme_classic(base_size = 20)


ggplot() +
  geom_density(
    aes(x = df_common$probability_with_mods.x, fill = "methylated transcripts"),
    alpha = 0.4,
    linewidth = 1
  ) +
  geom_density(
    aes(x = df_common$probability_without_mods.x, fill = "unmethylated transcripts"),
    alpha = 0.4,
    linewidth = 1
  ) +
  scale_fill_manual(
    name   = "",
    values = c(
      "methylated transcripts"   = "#1f77b4",
      "unmethylated transcripts" = "#ff7f0e"
    )
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  theme_classic(base_size = 20) +
  labs(
    x = "Probability",
    y = "Density",
    title = "Do transcripts end within 100nt of m6A sites? 28C1"
  )

ggplot() +
  geom_density(
    aes(x = df_common$probability_with_mods.y, fill = "methylated transcripts"),
    alpha = 0.4,
    linewidth = 1
  ) +
  geom_density(
    aes(x = df_common$probability_without_mods.y, fill = "unmethylated transcripts"),
    alpha = 0.4,
    linewidth = 1
  ) +
  scale_fill_manual(
    name   = "",
    values = c(
      "methylated transcripts"   = "#1f77b4",
      "unmethylated transcripts" = "#ff7f0e"
    )
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  theme_classic(base_size = 20) +
  labs(
    x = "Probability",
    y = "Density",
    title = "Do transcripts end within 100nt of m6A sites? 28K1"
  )

ggplot() +
  geom_density(
    aes(x = df_common$probability_with_mods.x, fill = "control"),
    alpha = 0.4,
    linewidth = 1
  ) +
  geom_density(
    aes(x = df_common$probability_with_mods.y, fill = "ks"),
    alpha = 0.4,
    linewidth = 1
  ) +
  scale_fill_manual(
    name   = "",
    values = c(
      "control"   = "#1f77b4",
      "ks" = "#ff7f0e"
    )
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  theme_classic(base_size = 20) +
  labs(
    x = "Probability",
    y = "Density",
    title = "Do transcripts end within 100nt of m6A sites?"
  )

ggplot() +
  geom_density(
    aes(x = df_common$probability_without_mods.x, fill = "control"),
    alpha = 0.4,
    linewidth = 1
  ) +
  geom_density(
    aes(x = df_common$probability_without_mods.y, fill = "ks"),
    alpha = 0.4,
    linewidth = 1
  ) +
  scale_fill_manual(
    name   = "",
    values = c(
      "control"   = "#1f77b4",
      "ks" = "#ff7f0e"
    )
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  theme_classic(base_size = 20) +
  labs(
    x = "Probability",
    y = "Density",
    title = "Do transcripts end within 100nt of m6A sites?"
  )

wilcox.test(df_common$probability_with_mods.x, df_common$probability_with_mods.y)
t.test(df_common$probability_with_mods.x, df_common$probability_with_mods.y)


# 
# 
# num_bins <- 30
# 
# p <- ggplot(df, aes(x = probability_with_mods,
#                     y = probability_without_mods)) +
#   geom_bin2d(bins = num_bins) +
#   geom_point(alpha = 0) +   # 👈 invisible points for ggMarginal
#   scale_x_continuous(
#     labels = percent_format(accuracy = 1)
#   ) +
#   scale_y_continuous(
#     labels = percent_format(accuracy = 1)
#   ) +
#   scale_fill_viridis_c(
#     trans = "log10",
#     name = "num sites"
#   ) +
#   labs(
#     y = "Probability (unmethylated)",
#     x = "Probability (methylated)"
#   ) +
#   guides(fill = "none") +
#   theme_classic(base_size = 20)
# 
# ggMarginal(
#   p,
#   type = "histogram",
#   bins = num_bins,
#   # fill = "grey70",
#   # color = "black"
# )

# ----- scatter: num_tx_within vs probability ----- #
# ggplot() +
#   geom_point(data = df, aes(x = probability_without_mods, y = percent_mod), alpha = 1, size=1) +
#   scale_y_continuous(
#     trans = logit_trans(),
#     breaks  = p_breaks,
#     labels  = percent_format(accuracy = 1),
#     name = "m6A / A (logit)"
#   ) +
#   scale_x_continuous(
#     trans = logit_trans(),
#     breaks  = p_breaks,
#     labels  = percent_format(accuracy = 1),
#     name = "Probability (logit)"
#   ) +
#   theme_classic(base_size = 20)

# ----- 2d heatmap: probability vs percent mod (logit) ----- #
# eps <- 1e-4   # or 1e-3 if you want tighter
# 
# df2 <- df %>%
#   mutate(
#     probability_without_mods = pmin(pmax(probability_without_mods, eps), 1 - eps),
#     probability_with_mods = pmin(pmax(probability_with_mods, eps), 1 - eps),
#     percent_mod = pmin(pmax(percent_mod, eps), 1 - eps)
#   )
# 
# ggplot(df2, aes(x = probability_with_mods, y = percent_mod)) +
#   geom_bin2d(bins = 30) +
#   scale_x_continuous(
#     trans = logit_trans(),
#     breaks  = p_breaks,
#     labels  = percent_format(accuracy = 1),
#     limits = c(eps, 1-eps),
#     name = "Probability (logit)"
#   ) +
#   scale_y_continuous(
#     trans = logit_trans(),
#     breaks  = p_breaks,
#     labels  = percent_format(accuracy = 1),
#     limits = c(min(df$percent_mod), 1 - eps),
#     name = "m6A / A"
#   ) +
#   scale_fill_viridis_c(
#     trans = "log2",
#     name = "num sites"
#   ) +
#   theme_classic(base_size = 20)

  