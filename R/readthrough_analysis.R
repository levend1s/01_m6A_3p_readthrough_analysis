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

cohensd_to_signif <- function(cd) {
  if (abs(cd) < 0.2) {
    return("ns")
  } else if (abs(cd) < 0.5) {
    return("*")
  } else if (abs(cd) < 0.8) {
    return("**")
  } else {
    return("***")
  }
}


library(ggplot2)
library(introdataviz)
library(ggsignif)
library(effectsize)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
tes_analysis_28hpi <- args[1]
tes_analysis_32hpi <- args[2]
tes_analysis_36hpi <- args[3]
out_file <- args[4]

ext <- file_ext(out_file)
name_no_ext <- file_path_sans_ext(out_file)

df_28 <- read.delim(tes_analysis_28hpi, header=TRUE, row.names="ID")
df_32 <- read.delim(tes_analysis_32hpi, header=TRUE, row.names="ID")
df_36 <- read.delim(tes_analysis_36hpi, header=TRUE, row.names="ID")


# -------------------- READTHROUGH proportion --------------------  #

group1 <- 'mean_rt_prop_g1'
group2 <- 'mean_rt_prop_g2'
y_label <- "run-on transcripts (%)"

x_label <- "hours post invasion (hpi)"
group_column_names <- c(group1, group2)

target_df_28 <- df_28[df_28$canonical_pa_site != 0, group_column_names]
target_df_32 <- df_32[df_32$canonical_pa_site != 0, group_column_names]
target_df_36 <- df_36[df_36$canonical_pa_site != 0, group_column_names]

# target_df_28 <- log(target_df_28)
# target_df_32 <- log(target_df_32)
# target_df_36 <- log(target_df_36)

tt_28 <- wilcox.test(target_df_28[[group1]], target_df_28[[group2]], alternative = "two.sided", paired=TRUE)
cd_28 <- cohens_d(target_df_28[[group1]], target_df_28[[group2]], paired = TRUE)
p_value_28 <- tt_28$p.value

tt_32 <- wilcox.test(target_df_32[[group1]], target_df_32[[group2]], alternative = "two.sided", paired=TRUE)
cd_32 <- cohens_d(target_df_32[[group1]], target_df_32[[group2]], paired = TRUE)
mean_change_32 <- mean(target_df_32[[group2]]) - mean(target_df_32[[group1]])
p_value_32 <- tt_32$p.value

tt_36 <- wilcox.test(target_df_36[[group1]], target_df_36[[group2]], alternative = "two.sided", paired=TRUE)
cd_36 <- cohens_d(target_df_36[[group1]], target_df_36[[group2]], paired = TRUE)
p_value_36 <- tt_36$p.value

#cohens_d(target_df_28[[group1]], target_df_28[[group2]], paired = TRUE)

num_points_28 <- nrow(target_df_28)
num_points_32 <- nrow(target_df_32)
num_points_36 <- nrow(target_df_36)

max_ys <- c(max(target_df_28) + max(target_df_28) * 0.1, max(target_df_32) + max(target_df_32) * 0.1, max(target_df_36) + max(target_df_36) * 0.1)
lims <- c(0, max(max_ys) * 1.1)

my_data = data.frame(
  y=c(target_df_28[,c(1)], target_df_28[,c(2)], target_df_32[,c(1)], target_df_32[,c(2)], target_df_36[,c(1)], target_df_36[,c(2)]),
  x=c(rep('28 HPI', num_points_28 * 2), rep('32 HPI', num_points_32 * 2), rep('36 HPI', num_points_36 * 2)),
  m=c(rep('control', num_points_28), rep('knock-sideways', num_points_28), rep('control', num_points_32), rep('knock-sideways', num_points_32), rep('control', num_points_36), rep('knock-sideways', num_points_36))
)

p <- ggplot(my_data, aes(x = x, y = y, fill = m)) +
  geom_split_violin() +
  ylim(lims) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outliers=FALSE) +
  stat_summary(fun.data = "mean_se", geom = "crossbar", show.legend = F, width=0.18, position = position_dodge(.2)) +
  geom_signif(
    y_position = max_ys, 
    xmin = c(0.75, 1.75, 2.75),
    xmax = c(1.25, 2.25, 3.25),
    # annotation = c(
    #   sprintf('p: %s\nc: %s', pval_to_signif(p_value_28), cohensd_to_signif(cd_28$Cohens_d)), 
    #   sprintf('p: %s\nc: %s', pval_to_signif(p_value_32), cohensd_to_signif(cd_32$Cohens_d)), 
    #   sprintf('p: %s\nc: %s', pval_to_signif(p_value_36), cohensd_to_signif(cd_36$Cohens_d))
    # ), 
    annotation = c(pval_to_signif(p_value_28), pval_to_signif(p_value_32), pval_to_signif(p_value_36)),
    tip_length = 0.005, 
    textsize = 10,
    size = 1
  ) +
  theme_classic(base_size = 20) +
  # theme(
  #   # legend.position = c(1, 1),                # top right (inside plot area)
  #   legend.justification = c(1, 0.5),           # anchor to top right
  #   legend.key.size = unit(1, "lines"),     # smaller legend keys
  #   legend.text = element_text(size = 16),     # smaller legend text
  #   legend.title = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   panel.grid.major = element_blank()
  # ) +
  # geom_line(aes(group = m), color = "gray70", alpha = 0.5) + # connect pairs
  scale_x_discrete(name = x_label, labels = c(sprintf('28 hpi\n(n=%d)', num_points_28), sprintf('32 hpi\n(n=%d)', num_points_32), sprintf('36 hpi\n(n=%d)', num_points_36))) +
  scale_y_continuous(name = y_label, limits=lims, expand = c(0, 0)) + 
  scale_fill_brewer(palette = "Set3", name = "") +
  theme(legend.position = "none")

ggsave(filename=paste0(name_no_ext, "_readthrough.", ext), plot = p, device=ext, dpi=300, width=6, height=4.8, units="in")







# -------------------- Poly_A --------------------  #


group1 <- 'mean_poly_a_length_g1'
group2 <- 'mean_poly_a_length_g2'
y_label <- "poly(A) length (nt)"

x_label <- "hours post invasion (hpi)"
group_column_names <- c(group1, group2)

target_df_28 <- df_28[df_28$canonical_pa_site != 0, group_column_names]
target_df_32 <- df_32[df_32$canonical_pa_site != 0, group_column_names]
target_df_36 <- df_36[df_36$canonical_pa_site != 0, group_column_names]

# target_df_28 <- log(target_df_28)
# target_df_32 <- log(target_df_32)
# target_df_36 <- log(target_df_36)

tt_28 <- wilcox.test(target_df_28[[group1]], target_df_28[[group2]], alternative = "two.sided", paired=TRUE)
cd_28 <- cohens_d(target_df_28[[group1]], target_df_28[[group2]], paired = TRUE)
p_value_28 <- tt_28$p.value

tt_32 <- wilcox.test(target_df_32[[group1]], target_df_32[[group2]], alternative = "two.sided", paired=TRUE)
cd_32 <- cohens_d(target_df_32[[group1]], target_df_32[[group2]], paired = TRUE)
mean_change_32 <- mean(target_df_32[[group2]]) - mean(target_df_32[[group1]])
p_value_32 <- tt_32$p.value

tt_36 <- wilcox.test(target_df_36[[group1]], target_df_36[[group2]], alternative = "two.sided", paired=TRUE)
cd_36 <- cohens_d(target_df_36[[group1]], target_df_36[[group2]], paired = TRUE)
p_value_36 <- tt_36$p.value

#cohens_d(target_df_28[[group1]], target_df_28[[group2]], paired = TRUE)

num_points_28 <- nrow(target_df_28)
num_points_32 <- nrow(target_df_32)
num_points_36 <- nrow(target_df_36)

max_ys <- c(max(target_df_28) + max(target_df_28) * 0.1, max(target_df_32) + max(target_df_32) * 0.1, max(target_df_36) + max(target_df_36) * 0.1)
lims <- c(0, max(max_ys) * 1.1)

my_data = data.frame(
  y=c(target_df_28[,c(1)], target_df_28[,c(2)], target_df_32[,c(1)], target_df_32[,c(2)], target_df_36[,c(1)], target_df_36[,c(2)]),
  x=c(rep('28 HPI', num_points_28 * 2), rep('32 HPI', num_points_32 * 2), rep('36 HPI', num_points_36 * 2)),
  m=c(rep('control', num_points_28), rep('knock-sideways', num_points_28), rep('control', num_points_32), rep('knock-sideways', num_points_32), rep('control', num_points_36), rep('knock-sideways', num_points_36))
)

p <- ggplot(my_data, aes(x = x, y = y, fill = m)) +
  geom_split_violin() +
  ylim(lims) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outliers=FALSE) +
  stat_summary(fun.data = "mean_se", geom = "crossbar", show.legend = F, width=0.18, position = position_dodge(.2)) +
  geom_signif(
    y_position = max_ys, 
    xmin = c(0.75, 1.75, 2.75),
    xmax = c(1.25, 2.25, 3.25),
    # annotation = c(
    #   sprintf('p: %s\nc: %s', pval_to_signif(p_value_28), cohensd_to_signif(cd_28$Cohens_d)), 
    #   sprintf('p: %s\nc: %s', pval_to_signif(p_value_32), cohensd_to_signif(cd_32$Cohens_d)), 
    #   sprintf('p: %s\nc: %s', pval_to_signif(p_value_36), cohensd_to_signif(cd_36$Cohens_d))
    # ), 
    annotation = c(pval_to_signif(p_value_28), pval_to_signif(p_value_32), pval_to_signif(p_value_36)),
    tip_length = 0.005, 
    textsize = 10,
    size = 1
  ) +
  theme_classic(base_size = 20) +
  # theme(
  #   # legend.position = c(1, 1),                # top right (inside plot area)
  #   legend.justification = c(1, 0.5),           # anchor to top right
  #   legend.key.size = unit(1, "lines"),     # smaller legend keys
  #   legend.text = element_text(size = 16),     # smaller legend text
  #   legend.title = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   panel.grid.major = element_blank()
  # ) +
  # geom_line(aes(group = m), color = "gray70", alpha = 0.5) + # connect pairs
  scale_x_discrete(name = x_label, labels = c(sprintf('28 hpi\n(n=%d)', num_points_28), sprintf('32 hpi\n(n=%d)', num_points_32), sprintf('36 hpi\n(n=%d)', num_points_36))) +
  scale_y_continuous(name = y_label, limits=lims, expand = c(0, 0)) + 
  scale_fill_brewer(palette = "Set3", name = "") +
  theme(legend.position = "none")

ggsave(filename=paste0(name_no_ext, "_poly_a_length.", ext), plot = p, device=ext, dpi=300, width=6, height=4.8, units="in")










# -------------------- UTR length --------------------  #


group1 <- 'mean_utr_length_g1'
group2 <- 'mean_utr_length_g2'
y_label <- "UTR length (bp)"

x_label <- "hours post invasion (hpi)"
group_column_names <- c(group1, group2)

target_df_28 <- df_28[df_28$canonical_pa_site != 0, group_column_names]
target_df_32 <- df_32[df_32$canonical_pa_site != 0, group_column_names]
target_df_36 <- df_36[df_36$canonical_pa_site != 0, group_column_names]

# target_df_28 <- log(target_df_28)
# target_df_32 <- log(target_df_32)
# target_df_36 <- log(target_df_36)

tt_28 <- wilcox.test(target_df_28[[group1]], target_df_28[[group2]], alternative = "two.sided", paired=TRUE)
cd_28 <- cohens_d(target_df_28[[group1]], target_df_28[[group2]], paired = TRUE)
p_value_28 <- tt_28$p.value

tt_32 <- wilcox.test(target_df_32[[group1]], target_df_32[[group2]], alternative = "two.sided", paired=TRUE)
cd_32 <- cohens_d(target_df_32[[group1]], target_df_32[[group2]], paired = TRUE)
mean_change_32 <- mean(target_df_32[[group2]]) - mean(target_df_32[[group1]])
p_value_32 <- tt_32$p.value

tt_36 <- wilcox.test(target_df_36[[group1]], target_df_36[[group2]], alternative = "two.sided", paired=TRUE)
cd_36 <- cohens_d(target_df_36[[group1]], target_df_36[[group2]], paired = TRUE)
p_value_36 <- tt_36$p.value

#cohens_d(target_df_28[[group1]], target_df_28[[group2]], paired = TRUE)

num_points_28 <- nrow(target_df_28)
num_points_32 <- nrow(target_df_32)
num_points_36 <- nrow(target_df_36)

max_ys <- c(max(target_df_28) + max(target_df_28) * 0.1, max(target_df_32) + max(target_df_32) * 0.1, max(target_df_36) + max(target_df_36) * 0.1)
lims <- c(0, max(max_ys) * 1.1)

my_data = data.frame(
  y=c(target_df_28[,c(1)], target_df_28[,c(2)], target_df_32[,c(1)], target_df_32[,c(2)], target_df_36[,c(1)], target_df_36[,c(2)]),
  x=c(rep('28 HPI', num_points_28 * 2), rep('32 HPI', num_points_32 * 2), rep('36 HPI', num_points_36 * 2)),
  m=c(rep('control', num_points_28), rep('knock-sideways', num_points_28), rep('control', num_points_32), rep('knock-sideways', num_points_32), rep('control', num_points_36), rep('knock-sideways', num_points_36))
)

p <- ggplot(my_data, aes(x = x, y = y, fill = m)) +
  geom_split_violin() +
  ylim(lims) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outliers=FALSE) +
  stat_summary(fun.data = "mean_se", geom = "crossbar", show.legend = F, width=0.18, position = position_dodge(.2)) +
  geom_signif(
    y_position = max_ys, 
    xmin = c(0.75, 1.75, 2.75),
    xmax = c(1.25, 2.25, 3.25),
    # annotation = c(
    #   sprintf('p: %s\nc: %s', pval_to_signif(p_value_28), cohensd_to_signif(cd_28$Cohens_d)), 
    #   sprintf('p: %s\nc: %s', pval_to_signif(p_value_32), cohensd_to_signif(cd_32$Cohens_d)), 
    #   sprintf('p: %s\nc: %s', pval_to_signif(p_value_36), cohensd_to_signif(cd_36$Cohens_d))
    # ), 
    annotation = c(pval_to_signif(p_value_28), pval_to_signif(p_value_32), pval_to_signif(p_value_36)),
    tip_length = 0.005, 
    textsize = 10,
    size = 1
  ) +
  theme_classic(base_size = 20) +
  # theme(
  #   # legend.position = c(1, 1),                # top right (inside plot area)
  #   legend.justification = c(1, 0.5),           # anchor to top right
  #   legend.key.size = unit(1, "lines"),     # smaller legend keys
  #   legend.text = element_text(size = 16),     # smaller legend text
  #   legend.title = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   panel.grid.major = element_blank()
  # ) +
  # geom_line(aes(group = m), color = "gray70", alpha = 0.5) + # connect pairs
  scale_x_discrete(name = x_label, labels = c(sprintf('28 hpi\n(n=%d)', num_points_28), sprintf('32 hpi\n(n=%d)', num_points_32), sprintf('36 hpi\n(n=%d)', num_points_36))) +
  scale_y_continuous(name = y_label, limits=lims, expand = c(0, 0)) + 
  scale_fill_brewer(palette = "Set3", name = "") +
  theme(legend.position = "none")

ggsave(filename=paste0(name_no_ext, "_utr_length.", ext), plot = p, device=ext, dpi=300, width=6, height=4.8, units="in")


