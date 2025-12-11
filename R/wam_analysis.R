# paired t test
# using for plot change in canonical, non canonical, and total wam
# also use for plotting proportion in num readthroughs
# these bars will also show the spread of the data
# maybe show them as violins? then we can see the underlying density

# install.packages("devtools")
# devtools::install_github("psyteachr/introdataviz")

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


library(ggplot2)
library(introdataviz)
library(ggsignif)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
wam_analysis_28hpi <- args[1]
wam_analysis_32hpi <- args[2]
wam_analysis_36hpi <- args[3]
out_file <- args[4]

ext <- file_ext(out_file)
name_no_ext <- file_path_sans_ext(out_file)

df_28 <- read.delim(wam_analysis_28hpi, header=TRUE, row.names="ID")
df_32 <- read.delim(wam_analysis_32hpi, header=TRUE, row.names="ID")
df_36 <- read.delim(wam_analysis_36hpi, header=TRUE, row.names="ID")

# -------------------- TOTAL M6A/A ratio --------------------  #
target_df_28 <- df_28[, c('total_wam_g1', 'total_wam_g2')] * 100
target_df_28 <- target_df_28[target_df_28$total_wam_g1 != 0, ]
target_df_28 <- target_df_28[target_df_28$total_wam_g2 != 0, ]

target_df_32 <- df_32[, c('total_wam_g1', 'total_wam_g2')] * 100
target_df_32 <- target_df_32[target_df_32$total_wam_g1 != 0, ]
target_df_32 <- target_df_32[target_df_32$total_wam_g2 != 0, ]

target_df_36 <- df_36[, c('total_wam_g1', 'total_wam_g2')] * 100
target_df_36 <- target_df_36[target_df_36$total_wam_g1 != 0, ]
target_df_36 <- target_df_36[target_df_36$total_wam_g2 != 0, ]
max_ys <- c(max(target_df_28) + max(target_df_28) * 0.1, max(target_df_32) + max(target_df_32) * 0.1, max(target_df_36) + max(target_df_36) * 0.1)
lims <- c(0, max(max_ys) * 1.1)
y_title <- expression(m^6 * "A / A (%)")
tt_28 <- t.test(target_df_28$total_wam_g1, target_df_28$total_wam_g2, paired = TRUE)
p_value_28 <- tt_28$p.value

tt_32 <- t.test(target_df_32$total_wam_g1, target_df_32$total_wam_g2, paired = TRUE)
p_value_32 <- tt_32$p.value

tt_36 <- t.test(target_df_36$total_wam_g1, target_df_36$total_wam_g2, paired = TRUE)
p_value_36 <- tt_36$p.value

num_points_28 <- nrow(target_df_28)
num_points_32 <- nrow(target_df_32)
num_points_36 <- nrow(target_df_36)

my_data = data.frame(
  y=c(target_df_28[,c(1)], target_df_28[,c(2)], target_df_32[,c(1)], target_df_32[,c(2)], target_df_36[,c(1)], target_df_36[,c(2)]),
  x=c(rep('28 HPI', num_points_28 * 2), rep('32 HPI', num_points_32 * 2), rep('36 HPI', num_points_36 * 2)),
  m=c(rep('control', num_points_28), rep('knock-sideways', num_points_28), rep('control', num_points_32), rep('knock-sideways', num_points_32), rep('control', num_points_36), rep('knock-sideways', num_points_36))
)

max_ys <- c(max(target_df_28), max(target_df_32), max(target_df_36))

# TODO go back to gene_methylation_analysis in rqc and figure out why some genes are showing up as >100% methylated!! fix that there not here!

total_plot <- ggplot(my_data, aes(x = x, y = y, fill = m)) +
  geom_split_violin() +
  # ylim(lims) +
  geom_boxplot(width = .2, fatten = TRUE, show.legend = FALSE, outliers = FALSE, coef=1) +
  # stat_summary(fun.data = "mean_se", geom = "crossbar", show.legend = F, width=0.18,
  #             position = position_dodge(.2)) +
  geom_signif(
    y_position = max_ys,
    xmin = c(0.75, 1.75, 2.75),
    xmax = c(1.25, 2.25, 3.25),
    annotation = c(pval_to_signif(p_value_28), pval_to_signif(p_value_32), pval_to_signif(p_value_36)),
    tip_length = 0.005,
    textsize = 7,
    size = 1
  ) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  # geom_line(aes(group = m), color = "gray70", alpha = 0.5) + # connect pairs
  scale_x_discrete(name = "", labels = c(sprintf('28 hpi\n(n=%d)', num_points_28), sprintf('32 hpi\n(n=%d)', num_points_32), sprintf('36 hpi\n(n=%d)', num_points_36))) +
  scale_y_continuous(name = y_title, limits = lims, expand = c(0, 0)) + 
  scale_fill_manual(values = c("#1bb6bb", "#f46b64"))

ggsave(filename=paste0(name_no_ext, "_total_wam.", ext), plot = total_plot, device=ext, dpi=300, width=6, height=4.8, units="in")






# -------------------- Canonical M6A/A ratio --------------------  #

# TODO go back to gene_methylation_analysis in rqc and figure out why some genes are showing up as >100% methylated!! fix that there not here!
target_df_28 <- df_28[, c('canonical_wam_g1', 'canonical_wam_g2')] * 100
target_df_28 <- target_df_28[
  target_df_28$canonical_wam_g1 > 0 & target_df_28$canonical_wam_g1 <= 100 &
  target_df_28$canonical_wam_g2 > 0 & target_df_28$canonical_wam_g2 <= 100,
]

target_df_32 <- df_32[, c('canonical_wam_g1', 'canonical_wam_g2')] * 100
target_df_32 <- target_df_32[
  target_df_32$canonical_wam_g1 > 0 & target_df_32$canonical_wam_g1 <= 100 &
    target_df_32$canonical_wam_g2 > 0 & target_df_32$canonical_wam_g2 <= 100,
]

target_df_36 <- df_36[, c('canonical_wam_g1', 'canonical_wam_g2')] * 100
target_df_36 <- target_df_36[
  target_df_36$canonical_wam_g1 > 0 & target_df_36$canonical_wam_g1 <= 100 &
    target_df_36$canonical_wam_g2 > 0 & target_df_36$canonical_wam_g2 <= 100,
]
lims <- c(0, 110)
max_ys <- c(100, 100, 100)
y_title <- expression(m^6 * "A / A (%)")
tt_28 <- t.test(target_df_28$canonical_wam_g1, target_df_28$canonical_wam_g2, paired = TRUE)
p_value_28 <- tt_28$p.value

tt_32 <- t.test(target_df_32$canonical_wam_g1, target_df_32$canonical_wam_g2, paired = TRUE)
p_value_32 <- tt_32$p.value

tt_36 <- t.test(target_df_36$canonical_wam_g1, target_df_36$canonical_wam_g2, paired = TRUE)
p_value_36 <- tt_36$p.value

num_points_28 <- nrow(target_df_28)
num_points_32 <- nrow(target_df_32)
num_points_36 <- nrow(target_df_36)

my_data = data.frame(
  y=c(target_df_28[,c(1)], target_df_28[,c(2)], target_df_32[,c(1)], target_df_32[,c(2)], target_df_36[,c(1)], target_df_36[,c(2)]),
  x=c(rep('28 HPI', num_points_28 * 2), rep('32 HPI', num_points_32 * 2), rep('36 HPI', num_points_36 * 2)),
  m=c(rep('control', num_points_28), rep('knock-sideways', num_points_28), rep('control', num_points_32), rep('knock-sideways', num_points_32), rep('control', num_points_36), rep('knock-sideways', num_points_36))
)

max_ys <- c(max(target_df_28), max(target_df_32), max(target_df_36))

canonical_plot <- ggplot(my_data, aes(x = x, y = y, fill = m)) +
  geom_split_violin() +
  # ylim(lims) +
  geom_boxplot(width = .2, fatten = TRUE, show.legend = FALSE, outliers = FALSE, coef=1) +
  # stat_summary(fun.data = "mean_se", geom = "crossbar", show.legend = F, width=0.18,
  #             position = position_dodge(.2)) +
  geom_signif(
    y_position = max_ys,
    xmin = c(0.75, 1.75, 2.75),
    xmax = c(1.25, 2.25, 3.25),
    annotation = c(pval_to_signif(p_value_28), pval_to_signif(p_value_32), pval_to_signif(p_value_36)),
    tip_length = 0.005,
    textsize = 7,
    size = 1
  ) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  # geom_line(aes(group = m), color = "gray70", alpha = 0.5) + # connect pairs
  scale_x_discrete(name = "", labels = c(sprintf('28 hpi\n(n=%d)', num_points_28), sprintf('32 hpi\n(n=%d)', num_points_32), sprintf('36 hpi\n(n=%d)', num_points_36))) +
  scale_y_continuous(name = y_title, limits = lims, expand = c(0, 0)) + 
  scale_fill_manual(values = c("#1bb6bb", "#f46b64"))

ggsave(filename=paste0(name_no_ext, "_canonical_wam.", ext), plot = canonical_plot, device=ext, dpi=300, width=6, height=4.8, units="in")









# -------------------- Non-Canonical M6A/A ratio --------------------  #
target_df_28 <- df_28[, c('non_canonical_wam_g1', 'non_canonical_wam_g2')] * 100
target_df_28 <- target_df_28[target_df_28$non_canonical_wam_g1 != 0, ]
target_df_28 <- target_df_28[target_df_28$non_canonical_wam_g2 != 0, ]

target_df_32 <- df_32[, c('non_canonical_wam_g1', 'non_canonical_wam_g2')] * 100
target_df_32 <- target_df_32[target_df_32$non_canonical_wam_g1 != 0, ]
target_df_32 <- target_df_32[target_df_32$non_canonical_wam_g2 != 0, ]

target_df_36 <- df_36[, c('non_canonical_wam_g1', 'non_canonical_wam_g2')] * 100
target_df_36 <- target_df_36[target_df_36$non_canonical_wam_g1 != 0, ]
target_df_36 <- target_df_36[target_df_36$non_canonical_wam_g2 != 0, ]

max_ys <- c(max(target_df_28) + max(target_df_28) * 0.1, max(target_df_32) + max(target_df_32) * 0.1, max(target_df_36) + max(target_df_36) * 0.1)
lims <- c(0, max(max_ys) * 1.1)
y_title <- expression(m^6 * "A / A (%)")
tt_28 <- t.test(target_df_28$non_canonical_wam_g1, target_df_28$non_canonical_wam_g2, paired = TRUE)
p_value_28 <- tt_28$p.value

tt_32 <- t.test(target_df_32$non_canonical_wam_g1, target_df_32$non_canonical_wam_g2, paired = TRUE)
p_value_32 <- tt_32$p.value

tt_36 <- t.test(target_df_36$non_canonical_wam_g1, target_df_36$non_canonical_wam_g2, paired = TRUE)
p_value_36 <- tt_36$p.value

num_points_28 <- nrow(target_df_28)
num_points_32 <- nrow(target_df_32)
num_points_36 <- nrow(target_df_36)

my_data = data.frame(
  y=c(target_df_28[,c(1)], target_df_28[,c(2)], target_df_32[,c(1)], target_df_32[,c(2)], target_df_36[,c(1)], target_df_36[,c(2)]),
  x=c(rep('28 HPI', num_points_28 * 2), rep('32 HPI', num_points_32 * 2), rep('36 HPI', num_points_36 * 2)),
  m=c(rep('control', num_points_28), rep('knock-sideways', num_points_28), rep('control', num_points_32), rep('knock-sideways', num_points_32), rep('control', num_points_36), rep('knock-sideways', num_points_36))
)

max_ys <- c(max(target_df_28), max(target_df_32), max(target_df_36))

# TODO go back to gene_methylation_analysis in rqc and figure out why some genes are showing up as >100% methylated!! fix that there not here!

non_canonical_plot <- ggplot(my_data, aes(x = x, y = y, fill = m)) +
  geom_split_violin() +
  # ylim(lims) +
  geom_boxplot(width = .2, fatten = TRUE, show.legend = FALSE, outliers = FALSE, coef=1) +
  # stat_summary(fun.data = "mean_se", geom = "crossbar", show.legend = F, width=0.18,
  #             position = position_dodge(.2)) +
  geom_signif(
    y_position = max_ys,
    xmin = c(0.75, 1.75, 2.75),
    xmax = c(1.25, 2.25, 3.25),
    annotation = c(pval_to_signif(p_value_28), pval_to_signif(p_value_32), pval_to_signif(p_value_36)),
    tip_length = 0.005,
    textsize = 7,
    size = 1
  ) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  # geom_line(aes(group = m), color = "gray70", alpha = 0.5) + # connect pairs
  scale_x_discrete(name = "", labels = c(sprintf('28 hpi\n(n=%d)', num_points_28), sprintf('32 hpi\n(n=%d)', num_points_32), sprintf('36 hpi\n(n=%d)', num_points_36))) +
  scale_y_continuous(name = y_title, limits = lims, expand = c(0, 0)) + 
  scale_fill_manual(values = c("#1bb6bb", "#f46b64"))

ggsave(filename=paste0(name_no_ext, "_non_canonical_wam.", ext), plot = non_canonical_plot, device=ext, dpi=300, width=6, height=4.8, units="in")
