# show differential abundance of: 
# 1. expression counts
# 2. count of transcripts which overlap multiple genes (ambiguous reads)

# Show as:
# 1. volcano plot
# 2. MD plot

library(edgeR)
library(ggplot2)
library(pheatmap)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
out_file <- args[2]
columns_arg <- as.numeric(strsplit(args[3], ",")[[1]])

# --- normal DE (Fig 2C)
count_matrix <- input_file
x <- read.delim(count_matrix, header=TRUE, skip=1, row.names="Geneid")
columns <- x[,columns_arg] # 36hpi
group <- factor(c("control","control","knock-sideways","knock-sideways"))

# create matrix, filter genes with less than 10 reads, 
# normalize counts by library sizes, estimate dispersion across 
# replicates (ie handles variability between replicates)
y <- DGEList(counts=columns, group=group)
y <- y[filterByExpr(y, group=group, min.count=10), , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)

# Negative binomial GLM (edge v2)
fit <- glmFit(y, design)
qlf_nbglm <- glmLRT(fit, coef=2)
qlf <- qlf_nbglm

# plot!
tt = topTags(qlf,n=Inf)
df <- as.data.frame(tt)
colnames <- c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")

point_size <- 2

df$Significant <- "none"
df$Significant[df$logFC > 0 & df$FDR < 0.05] <- "up"
df$Significant[df$logFC < 0 & df$FDR < 0.05] <- "down"

df_not <- df[df$Significant == "none", ]
df_sig <- df[df$Significant != "none", ]

p <- ggplot() +
  geom_point(data = df_not, aes(x = logCPM, y = logFC), color = "grey", alpha = 1, size = point_size) +
  geom_point(data = df_sig, aes(x = logCPM, y = logFC, color = Significant), alpha = 1, size = point_size) +
  # geom_hline(yintercept = -log10(0.05), color = "red") +
  scale_color_manual(values = c("down" = "blue", "up" = "red", "none" = "gray")) +
  geom_point(alpha = 1.0) +
  theme_classic(base_size = 20) +
  labs(x = "logCPM", y = "logFC") +
  theme(legend.position = "none")

ext <- file_ext(out_file)

ggsave(
    filename = out_file,   # output filename
    plot = p,                  # ggplot object
    device = ext,            # EPS format
    width = 6,                 # width in inches
    height = 4.8,                # height in inches
    units = "in",              # units (inches, cm, mm)
    dpi = 300                  # resolution, optional (ignored by EPS)
)
