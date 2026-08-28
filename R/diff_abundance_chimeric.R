# show differential abundance of: 
# 1. expression counts
# 2. count of transcripts which overlap multiple genes (ambiguous reads)

# Show as:
# 1. volcano plot
# 2. MD plot

library(edgeR)
library(ggplot2)
library(tools)
library(GO.db)
library(clusterProfiler)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
chimeric_matrix <- args[1]
chimeric_columns <- as.numeric(strsplit(args[2], ",")[[1]])
transcript_matrix <- args[3]
transcript_columns <- as.numeric(strsplit(args[4], ",")[[1]])
gaf_file <- args[5]
out_file <- args[6]

# ---------- chimeric DE (interaction term between treatment AND chimeric/bulk count, figure 5C)
x_chimeric <- read.delim(chimeric_matrix, header=TRUE, row.names="Geneid")
x <- read.delim(transcript_matrix, header=TRUE, skip=1, row.names="Geneid")

x_chimeric <- x_chimeric[ !grepl("API|MIT", rownames(x_chimeric), ignore.case = TRUE), , drop = FALSE ]
x <- x[ !grepl("API|MIT", rownames(x), ignore.case = TRUE), , drop = FALSE ]

identical(rownames(x_chimeric), rownames(x))

columns_chimeric <- x_chimeric[,chimeric_columns]
columns <- x[,transcript_columns]

not_chimeric <- columns - columns_chimeric


merged_df <- cbind(columns_chimeric, not_chimeric)

group <- factor(rep(c("control", "knocksideways"), each=2, times=2))
chimeric_status <- factor(rep(c("chimeric", "not_chimeric"), each=4))

# explicit zero intercept
design <- model.matrix(~0 + group * chimeric_status)
y <- DGEList(counts=merged_df, group=group)

HasCoverage <- rowSums(columns >= 10) == 4
HasChimeric <- rowSums(columns_chimeric >= 10) >= 2

y <- y[HasCoverage & HasChimeric,, keep.lib.sizes=FALSE]

TotalLibSize <- y$samples$lib.size[chimeric_status=="chimeric"] + y$samples$lib.size[chimeric_status=="not_chimeric"]

y$samples$lib.size <- rep(TotalLibSize, times=2)

y <- estimateDisp(y, design=design, trend="none")
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = "groupknocksideways:chimeric_statusnot_chimeric")

tt = topTags(lrt,n=Inf)
df = as.data.frame(tt)

point_size <- 2
pval_cutoff <- 0.01

df$negLogFDR <- -log10(df$FDR)
df$logFC <- -df$logFC

df$Significant <- "none"
df$Significant[df$logFC > 0 & df$FDR < pval_cutoff] <- "up"
df$Significant[df$logFC < 0 & df$FDR < pval_cutoff] <- "down"

df_not <- df[df$Significant == "none", ]
df_sig <- df[df$Significant != "none", ]


# ---------- raw normal chimeric DE (no interaction term)
# group <- factor(c("control","control","knock-sideways","knock-sideways"))
# 
# # create matrix, filter genes with less than 10 reads,
# # normalize counts by library sizes, estimate dispersion across
# # replicates (ie handles variability between replicates)
# y <- DGEList(counts=columns_chimeric, group=group)
# # y <- y[filterByExpr(y, group=group, min.count=10), , keep.lib.sizes=FALSE]
# 
# y$samples$lib.size <- weights
# y <- calcNormFactors(y)
# 
# 
# design <- model.matrix(~group)
# y <- estimateDisp(y, design)
# 
# # Negative binomial GLM (edge v2)
# fit <- glmFit(y, design)
# qlf_nbglm <- glmLRT(fit, coef=2)
# qlf <- qlf_nbglm
# 
# # qlf <- exactTest(y)
# 
# # plot!
# tt = topTags(qlf,n=Inf)
# df <- as.data.frame(tt)
# colnames <- c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")
# 
# point_size <- 2
# 
# df$Significant <- "none"
# df$Significant[df$logFC > 0 & df$FDR < 0.05] <- "up"
# df$Significant[df$logFC < 0 & df$FDR < 0.05] <- "down"
# 
# df_not <- df[df$Significant == "none", ]
# df_sig <- df[df$Significant != "none", ]
# 
# p <- ggplot() +
#   geom_point(data = df_not, aes(x = logCPM, y = logFC), color = "grey", alpha = 1, size = point_size) +
#   geom_point(data = df_sig, aes(x = logCPM, y = logFC, color = Significant), alpha = 1, size = point_size) +
#   # geom_hline(yintercept = -log10(0.05), color = "red") +
#   scale_color_manual(values = c("down" = "blue", "up" = "red", "none" = "gray")) +
#   geom_point(alpha = 1.0) +
#   theme_classic(base_size = 20) +
#   labs(x = "logCPM", y = "logFC") +
#   geom_hline(yintercept = 0, color = "red") +
#   theme(legend.position = "none")
# p


# ----------  % chimeric before vs % chimeric after KS
percent_chimeric <- columns_chimeric / columns
percent_chimeric[is.na(percent_chimeric)] <- 0

colnames(columns) <- colnames(columns_chimeric)
weights <- colSums(columns)
labels <- colnames(columns)
avg_C <- rowSums(percent_chimeric[, c(labels[1], labels[2])] * weights[c(labels[1], labels[2])]) / sum(weights[c(labels[1], labels[2])])
avg_K <- rowSums(percent_chimeric[, c(labels[3], labels[4])] * weights[c(labels[3], labels[4])]) / sum(weights[c(labels[3], labels[4])])

avg_perc <- rowSums(percent_chimeric[, labels] * weights) / sum(weights)

avg_count_C <- rowSums(columns_chimeric[, c(labels[1], labels[2])] * weights[c(labels[1], labels[2])]) / sum(weights[c(labels[1], labels[2])])
avg_count_K <- rowSums(columns_chimeric[, c(labels[3], labels[4])] * weights[c(labels[3], labels[4])]) / sum(weights[c(labels[3], labels[4])])
result_percent <- data.frame(avg_C = avg_C, avg_K = avg_K, avg_count_C = avg_count_C, avg_count_K = avg_count_K, avg_perc=avg_perc)

min_val <- 1e-3
max_val <- 1 - min_val

cols <- c("avg_C", "avg_K", "avg_perc")#  "avg_count_C", "avg_count_K")  # columns to modify

result_percent[cols][result_percent[cols] <= min_val] <- min_val
result_percent[cols][result_percent[cols] >= 1-min_val] <- max_val

result_percent$logit_avg_C <- log(result_percent$avg_C / (1 - result_percent$avg_C))
result_percent$logit_avg_K <- log(result_percent$avg_K / (1 - result_percent$avg_K))
result_percent$logit_avg_perc <- log(result_percent$avg_perc / (1 - result_percent$avg_perc))

result_percent$percent_change <- result_percent$avg_K - result_percent$avg_C

# signif_df <- data.frame(Significant = df$Significant, row.names = rownames(df))
# Merge by row names
result_percent <- merge(result_percent, df, by = "row.names", all.y = TRUE)
rownames(result_percent) <- result_percent$Row.names
result_percent$Row.names <- NULL

breaks <- c(min_val, 0.01, 0.5, 0.99, max_val)
breaks_labels <- breaks
breaks_labels[c(1, length(breaks_labels))] <- c(0, 1)
log_breaks <- log(breaks / (1 - breaks))

result_percent$avg_count_C <- result_percent$avg_count_C + 1
result_percent$avg_count_K <- result_percent$avg_count_K + 1
# result_percent <- result_percent[HasCoverage & HasChimeric,]

df_not_result_percent <- result_percent[result_percent$Significant == "none", ]
df_sig_result_percent <- result_percent[result_percent$Significant != "none", ]

# ------------ write files ! ------------ #

ext <- file_ext(out_file)
name_no_ext <- file_path_sans_ext(out_file)

write.table(tt, file=paste0(name_no_ext, "_table.txt"))

# md
p <- ggplot() +
  geom_point(data = df_not_result_percent, aes(x = logit_avg_C, y = logFC), color = "grey") +
  geom_point(data = df_sig_result_percent, aes(x = logit_avg_C, y = logFC, color = Significant)) +
  # geom_hline(yintercept = -log10(0.05), color = "red") +
  scale_color_manual(values = c("down" = "blue", "up" = "red", "none" = "gray")) +
  geom_point(alpha = 1.0) +
  theme_classic(base_size = 20) +
  geom_hline(yintercept = 0, color = "red") +
  labs(x = "control chimeric / total (%)", y = "logFC") +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = log(breaks /
                   (1 - breaks)),
    labels = paste0(breaks_labels * 100, "%")
  )
p

md_file <- paste0(name_no_ext, "_md.", ext)

p <- p +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.title = element_blank()
  )

ggsave(
    filename = md_file,   # output filename
    plot = p,                  # ggplot object
    device = ext,            # EPS format
    width = 6,                 # width in inches
    height = 4.8,                # height in inches
    units = "in",              # units (inches, cm, mm)
    bg = "transparent",
    dpi = 300                  # resolution, optional (ignored by EPS)
)

# volcano
p <- ggplot() +
  geom_point(data = df_not, aes(x = logFC, y = negLogFDR), color = "grey", alpha = 1, size = point_size) +
  geom_point(data = df_sig, aes(x = logFC, y = negLogFDR, color = Significant), alpha = 1, size = point_size) +
  geom_hline(yintercept = -log10(pval_cutoff), color = "red") +
  scale_color_manual(values = c("down" = "blue", "up" = "red", "none" = "gray")) +
  geom_point(alpha = 1.0) +
  theme_classic(base_size = 20) +
  labs(x = "logFC", y = "-logFDR") +
  theme(legend.position = "none")

volcano_file <- paste0(name_no_ext, "_volcano.", ext)

ggsave(
    filename = volcano_file,   # output filename
    plot = p,                  # ggplot object
    device = ext,            # EPS format
    width = 6,                 # width in inches
    height = 4.8,                # height in inches
    units = "in",              # units (inches, cm, mm)
    dpi = 300                  # resolution, optional (ignored by EPS)
)


# perc vs perc
perc_plot <- ggplot() +
  # geom_point(alpha = 1.0) +
  geom_point(data = df_not_result_percent, aes(x = logit_avg_C, y = logit_avg_K), color = "grey", alpha = 1) +
  geom_point(data = df_sig_result_percent, aes(x = logit_avg_C, y = logit_avg_K, color = Significant), alpha = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_x_continuous(
    breaks = log_breaks,
    labels = paste0(breaks_labels * 100, "%"),
    limits = c(log_breaks[1], log_breaks[length(log_breaks)])
  ) +
  scale_y_continuous(
    breaks = log_breaks,
    labels = paste0(breaks_labels * 100, "%"),
    limits = c(log_breaks[1], log_breaks[length(log_breaks)])
  ) +
  scale_color_manual(values = c("down" = "blue", "up" = "red", "none" = "gray")) +
  labs(x = "% chimeric transcripts", y = "% chimeric transcripts") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

perc_file <- paste0(name_no_ext, "_perc_vs_perc.", ext)

ggsave(
  filename = perc_file,   # output filename
  plot = perc_plot,                  # ggplot object
  device = ext,            # EPS format
  width = 6,                 # width in inches
  height = 4.8,                # height in inches
  units = "in",              # units (inches, cm, mm)
  dpi = 300                  # resolution, optional (ignored by EPS)
)

# count plot
count_plot <- ggplot() +
  geom_point(data = df_not_result_percent, aes(x = avg_count_C, y = avg_count_K), color = "grey", alpha = 1) +
  geom_point(data = df_sig_result_percent, aes(x = avg_count_C, y = avg_count_K, color = Significant), alpha = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_x_continuous(
    limits = c(1, 10000),
    trans="log10"
  ) +
  scale_y_continuous(
    limits = c(1, 10000),
    trans="log10"
  ) +
  scale_color_manual(values = c("down" = "blue", "up" = "red", "none" = "gray")) +
  labs(x = "chimeric transcript count (log10(count + 1))", y = "chimeric transcript count (log10(count + 1))") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

count_file <- paste0(name_no_ext, "_count_vs_count.", ext)

ggsave(
  filename = count_file,   # output filename
  plot = count_plot,                  # ggplot object
  device = ext,            # EPS format
  width = 6,                 # width in inches
  height = 4.8,                # height in inches
  units = "in",              # units (inches, cm, mm)
  dpi = 300                  # resolution, optional (ignored by EPS)
)


# # GSEA

# # ---- GO ENRICHMENT
gaf <- read.delim(gaf_file, header = FALSE, comment.char = "!", stringsAsFactors = FALSE, sep="\t")
colnames(gaf)[c(2,5,10)] <- c("GeneID", "GO_ID", "Description")
term2gene <- gaf[, c("GO_ID", "GeneID")]
go_ids <- unique(term2gene$GO_ID)
go_names <- Term(go_ids)
term2name <- data.frame(GO_ID = go_ids, GO_Name = go_names)


# gois <- rownames(df[df$Significant!="none",])
# go_enrich <- enricher(
#   gene = gois,
#   TERM2GENE = term2gene,
#   TERM2NAME = term2name
# )
# dotplot(go_enrich)

gene_list <- sign(df$logFC) * df$negLogFDR
names(gene_list) <- rownames(df)

# Sort descending for GSEA
gene_list <- sort(gene_list, decreasing = TRUE)

# rank based list
go_gsea <- GSEA(
  geneList   = gene_list,
  TERM2GENE  = term2gene,
  TERM2NAME  = term2name,
  minGSSize  = 5,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# ------------ write files ! ------------ #

# dotplot(go_gsea)
num_enriched_terms <- nrow(go_gsea)

if (is.null(go_gsea) || nrow(go_gsea@result) == 0) {

    message("No enriched GO terms found.")

    # Save an empty placeholder file so downstream steps don't break
    empty_file <- paste0(name_no_ext, "_NO_ENRICHED_TERMS.", ext)

    # Simple blank plot or text message plot
    dotplot <- ggplot() + 
         theme_void() +
         annotate("text", x = 0.5, y = 0.5, label = "No enriched GO terms", size = 6)

    ridgeplot <- ggplot() + 
        theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No enriched GO terms", size = 6)
} else {
    dotplot <- dotplot(
        go_gsea,
        x = "NES",        # <- x-axis is number of genes
        size = "GeneRatio", # dot size = fraction of pathway genes in your list
        showCategory = num_enriched_terms,  # top 20 terms
        color = "p.adjust"  # color by adjusted p-value
    )

    ridgeplot <- ridgeplot(go_gsea)
}

ggsave(
    filename = paste0(name_no_ext, "_dotplot.", ext),   # output filename
    plot = dotplot,                  # ggplot object
    device = ext,            # EPS format
    # width = 6,                 # width in inches
    # height = 4.8,                # height in inches
    units = "in",              # units (inches, cm, mm)
    dpi = 300                  # resolution, optional (ignored by EPS)
)

ggsave(
    filename = paste0(name_no_ext, "_ridgeplot.", ext),   # output filename
    plot = ridgeplot,                  # ggplot object
    device = ext,            # EPS format
    # width = 6,                 # width in inches
    # height = 4.8,                # height in inches
    units = "in",              # units (inches, cm, mm)
    dpi = 300                  # resolution, optional (ignored by EPS)
)


