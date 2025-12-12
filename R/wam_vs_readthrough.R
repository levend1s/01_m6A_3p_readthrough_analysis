# tes vs wam analysis


library(ggplot2)
library(tools)


args <- commandArgs(trailingOnly = TRUE)
num_canonical_mods <- args[1]
gene_neighbour_list <- args[2]

tes_analysis_28hpi <- args[3]
tes_analysis_32hpi <- args[4]
tes_analysis_36hpi <- args[5]

wam_analysis_28hpi <- args[6]
wam_analysis_32hpi <- args[7]
wam_analysis_36hpi <- args[8]

out_file <- args[9]

ext <- file_ext(out_file)
name_no_ext <- file_path_sans_ext(out_file)



df_28_wam <- read.delim(wam_analysis_28hpi, header=TRUE)
df_32_wam <- read.delim(wam_analysis_32hpi, header=TRUE)
df_36_wam <- read.delim(wam_analysis_36hpi, header=TRUE)

df_28_tes <- read.delim(tes_analysis_28hpi, header=TRUE)
df_32_tes <- read.delim(tes_analysis_32hpi, header=TRUE)
df_36_tes <- read.delim(tes_analysis_36hpi, header=TRUE)

df_neighbours <- read.delim(gene_neighbour_list, header=FALSE)
colnames(df_neighbours) <- c("contig", "start", "end", "name", "score", "strand", "type")
df_neighbours <- df_neighbours[df_neighbours$type == "protein_coding_gene",]
df_neighbours$ID <- sub(".*ID=([^;]+);.*", "\\1", df_neighbours$name)

df_28 <- merge(df_28_wam, df_28_tes, by = "ID")
df_28$ID_clean <- sub("\\..*$", "", df_28$ID)
df_28$wam_diff <- df_28$canonical_wam_g1 - df_28$canonical_wam_g2
df_28$tes_diff <- df_28$mean_rt_prop_g2 - df_28$mean_rt_prop_g1
df_28$canonical_mods <- gsub("\\[|\\]", "", df_28$canonical_mods)
df_28$canonical_mods <- strsplit(df_28$canonical_mods, ",\\s*")
df_28$num_canonical_mods <- lengths(df_28$canonical_mods)

df_28 <- df_28[df_28$wam_diff != 0,]
df_28 <- df_28[df_28$tes_diff != 0,]

if (num_canonical_mods != "all") {
    df_28 <- df_28[df_28$num_canonical_mods == num_canonical_mods,]
}


df_28 <- df_28[df_28$ID_clean %in% df_neighbours$ID, ]
print(nrow(df_28))

target_df_28 <- df_28[df_28$canonical_pa_site != 0, c('wam_diff', 'tes_diff')] * 100

cor_val <- cor(target_df_28$wam_diff, target_df_28$tes_diff, method = "pearson")

nrows <- nrow(target_df_28)

p <- ggplot(target_df_28, aes(x = wam_diff, y = tes_diff)) +
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
  scale_y_continuous(name="change in run-ons (%)", trans = "identity") +
  scale_x_continuous(name=sprintf("change in methylation (%%)\n(n=%d)", nrows), expand = c(0, 0)) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = Inf, y = Inf,
           label = paste("r =", round(cor_val, 2)),
           hjust = 1, vjust = 2, size = 10)

ggsave(filename=paste0(name_no_ext, "_wam_vs_readthrough_correlation.", ext), plot = p, device=ext, dpi=300, width=6, height=4.8, units="in")
