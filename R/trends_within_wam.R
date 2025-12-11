# trends within WAM

# TODO plot hidden trends in my gene_methylation_analysis, eg % methylation vs read depth, % methylation vs num_canonical_mods, num_canonical_mods vs read depth
# I would think that genes with more canonical m6A sites are probably those higher expressed, and these have lower WAMs than genes with a single cm6A, since it appears
# some m6As are more methylated than others

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

# depth vs num canonical mods (expected positive)
# depth vs wam (expected nothing? or negative)
# wam vs num canonical mods (expected negative)

# conclusion, more depth means more poly adenylation sites
# more poly adenylation sites means more canonical mods
# more canonical mods means reduced WAM poser


# ------------------- depth vs wam ------------------- #
# convert canonical mods into list and count into num_canonical_mods
df_28$canonical_mods <- gsub("\\[|\\]", "", df_28$canonical_mods)
df_28$canonical_mods <- strsplit(df_28$canonical_mods, ",\\s*")
df_28$num_canonical_mods <- lengths(df_28$canonical_mods)

target_df_28 <- df_28
target_df_28$canonical_wam_g1 <-target_df_28$canonical_wam_g1 * 100
target_df_28$canonical_wam_g2 <-target_df_28$canonical_wam_g2 * 100
target_df_28 <- target_df_28[
  target_df_28$canonical_wam_g1 > 0 & target_df_28$canonical_wam_g1 <= 100 &
    target_df_28$canonical_wam_g2 > 0 & target_df_28$canonical_wam_g2 <= 100,
]

cor_val <- cor(target_df_28$canonical_wam_g1, target_df_28$average_depth_g1, method = "pearson")
# 
p <- ggplot(target_df_28, aes(x = canonical_wam_g1, y = average_depth_g1)) +
  geom_point(size=1) +
  theme_classic(base_size = 20) +
  theme(
    legend.justification = c(1, 0.5),           # anchor to top right
    legend.key.size = unit(1, "lines"),     # smaller legend keys
    legend.text = element_text(size = 16),     # smaller legend text
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  scale_y_continuous(name="max read depth (count)", trans = "log2") +
  scale_x_continuous(name="methylation (%)", expand = c(0, 0)) +
  # geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = Inf, y = Inf,
           label = paste("r =", round(cor_val, 2)),
           hjust = 1, vjust = 2, size = 10)

ggsave(filename=paste0(name_no_ext, "_depth_vs_wam.", ext), plot = p, device=ext, dpi=300, width=6, height=4.8, units="in")







# ------------------- depth vs wam ------------------- #

num_points <- nrow(target_df_28)

my_data = data.frame(
  y=c(target_df_28$average_depth_g1, target_df_28$average_depth_g2),
  x=as.character(c(target_df_28$num_canonical_mods, target_df_28$num_canonical_mods)),
  m=c(rep('control', num_points), rep('knock-sideways', num_points))
)
y_trans <- "log2"
y_label <- "max read depth (count)"

x_ticks <- seq(0, max(target_df_28$num_canonical_mods))

sprintf('28 hpi\n(n=%d)', num_points)
  
  
max_n <- max(target_df_28$num_canonical_mods)
min_n <- min(target_df_28$num_canonical_mods)
counts <- tabulate(target_df_28$num_canonical_mods, nbins = max_n)
counts <- as.integer(counts)
# counts[i] = number of times (i-1) occurs in x
x_axes_labels <- sprintf('%d\n(n=%d)', seq(min_n, max_n), counts)

p <- ggplot(my_data, aes(x = x, y = y, fill = m)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  theme(
    legend.justification = c(1, 0.5),           # anchor to top right
    legend.key.size = unit(1, "lines"),     # smaller legend keys
    legend.text = element_text(size = 16),     # smaller legend text
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  scale_y_continuous(name=y_label, trans = y_trans) +
  scale_x_discrete(name="canonical m6A (count)", labels = x_axes_labels)

ggsave(filename=paste0(name_no_ext, "_depth_vs_num_canonical_m6A.", ext), plot = p, device=ext, dpi=300, width=12, height=9.6, units="in")







# ------------------- depth vs wam ------------------- #

num_points <- nrow(target_df_28)

my_data = data.frame(
  y=c(target_df_28$canonical_wam_g1, target_df_28$canonical_wam_g2),
  x=as.character(c(target_df_28$num_canonical_mods, target_df_28$num_canonical_mods)),
  m=c(rep('control', num_points), rep('knock-sideways', num_points))
)
y_trans <- "identity"
y_label <- "methylation (%)"

x_ticks <- seq(0, max(target_df_28$num_canonical_mods))

sprintf('28 hpi\n(n=%d)', num_points)
  
  
max_n <- max(target_df_28$num_canonical_mods)
min_n <- min(target_df_28$num_canonical_mods)
counts <- tabulate(target_df_28$num_canonical_mods, nbins = max_n)
counts <- as.integer(counts)
# counts[i] = number of times (i-1) occurs in x
x_axes_labels <- sprintf('%d\n(n=%d)', seq(min_n, max_n), counts)

p <- ggplot(my_data, aes(x = x, y = y, fill = m)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  theme(
    legend.justification = c(1, 0.5),           # anchor to top right
    legend.key.size = unit(1, "lines"),     # smaller legend keys
    legend.text = element_text(size = 16),     # smaller legend text
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  scale_y_continuous(name=y_label, trans = y_trans) +
  scale_x_discrete(name="canonical m6A (count)", labels = x_axes_labels)

ggsave(filename=paste0(name_no_ext, "_wam_vs_num_canonical_m6A.", ext), plot = p, device=ext, dpi=300, width=12, height=9.6, units="in")
