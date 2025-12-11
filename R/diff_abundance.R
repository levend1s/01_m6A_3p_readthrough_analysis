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

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
columns_arg <- as.numeric(strsplit(args[2], ",")[[1]])
gaf_file <- args[3]
out_file <- args[4]


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
name_no_ext <- file_path_sans_ext(out_file)

md_file <- paste0(name_no_ext, "_md.", ext)

ggsave(
    filename = md_file,   # output filename
    plot = p,                  # ggplot object
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


# gois <- rownames(df[df$Significant=="none",])
# go_enrich <- enricher(
#   gene = gois,
#   TERM2GENE = term2gene,
#   TERM2NAME = term2name
# )
# dotplot(go_enrich)

head(df)
gene_list <- sign(df$logFC) * df$LR
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
# cnetplot(go_gsea)
# gseaplot(gsea_res, geneSetID=1)

