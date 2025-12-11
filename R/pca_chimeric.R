# show differential abundance of: 
# 1. expression counts
# 2. count of transcripts which overlap multiple genes (ambiguous reads)

# Show as:
# 1. volcano plot
# 2. MD plot

# BiocManager::install("clusterProfiler")
# BiocManager::install("AnnotationHub")
# BiocManager::install("KEGGREST")

library(edgeR)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(dplyr)
library(compositions)
library(KEGGREST)
library(GO.db)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
count_matrix <- args[1]
gaf_file <- args[2]
out_file <- args[3]

ext <- file_ext(out_file)
name_no_ext <- file_path_sans_ext(out_file)

# ------------- transcript count PCA (Fig 2A,B)
# 

x <- read.delim(count_matrix, header=TRUE, row.names="Geneid")
row_names <- c("28C1", "28C2", "28K1", "28K2", "32C1", "32C2", "32K1", "32K2", "36C1", "36C2", "36K1", "36K2")
colnames(x) <- row_names
columns <- x
group <- factor(rep(c("control","knocksideways"), each=2, times=3))
time  <- factor(rep(c("28","32","36"), each=4))

# create matrix, filter genes with less than 10 reads,
# normalize counts by library sizes, estimate dispersion across
# replicates (ie handles variability between replicates)
y <- DGEList(counts=columns, group=group)
y <- y[filterByExpr(y, group=group, min.count=10), , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group*time)
y <- estimateDisp(y, design)

# Negative binomial GLM (edge v2)
fit <- glmFit(y, design)
lrt1 <- glmLRT(fit, coef="groupknocksideways")

# does the treatment change differently over time?
lrt2 <- glmLRT(fit, coef=grep("groupknocksideways:time", colnames(design)))

# plot!
ttlrt1 = topTags(lrt1,n=Inf)
ttlrt2 = topTags(lrt2,n=Inf)
colnames <- c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")
plotMD(lrt1)

# logCPM from edgeR object
logCPM <- cpm(y, log=TRUE, prior.count=1)

# PCA
pca <- prcomp(t(logCPM), center = TRUE)
data <- logCPM
# num_clusters <- 6
row_names <- colnames(logCPM)
method <- "chimeric count"











# ---------- PLOT PCA 
pca_df <- data.frame(pca$x, group=group)
# assuming 'pca' is result from prcomp(t(norm_dat), scale.=TRUE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
var_exp <- round(var_exp[1:10] * 100, 2)     # % variance for first 10 PCs

# Scree plot
scree_df <- data.frame(PC = seq_along(var_exp), VarExp = var_exp)

screeplot <- ggplot(scree_df, aes(x = PC, y = VarExp)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "PC", y = "% Variance Explained")

ggsave(filename=paste0(name_no_ext, "_screeplot.", ext), plot = screeplot, device=ext, dpi=300)


pca_plot <- ggplot(pca_df, aes(x=PC1, y=PC2, color=group, shape=time)) +
  geom_point(size=4) +
  # stat_ellipse(aes(group = Cluster), type = "norm", level = 0.8) +
  labs(
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "right",     # legend outside
    legend.box = "vertical",
    plot.title = element_text(hjust=0.5)
  )

# export 740x480
ggsave(filename=paste0(name_no_ext, "_pca_plot.", ext), plot = pca_plot, device=ext, dpi=300, width=7.4, height=4.8, units="in")

# got top contributing genes IDs
loadings <- pca$rotation  # genes × PCs
num_top_genes <- 50#nrow(loadings)
top_PC1_genes <- names(sort(abs(loadings[,1]), decreasing=TRUE)[1:num_top_genes])
bottom_PC1_genes <- names(sort(abs(loadings[,1]), decreasing=FALSE)[1:num_top_genes])

all_pc1_genes <- sort(abs(loadings[,1]), decreasing=TRUE)
all_pc2_genes <- sort(abs(loadings[,2]), decreasing=TRUE)

names(all_pc1_genes) <- sub("\\.[0-9]+$", "", names(all_pc1_genes))
all_pc1_genes <- all_pc1_genes[!duplicated(names(all_pc1_genes))]

# Don't really care about PC2 cause it's not the main separator
top_PC2_genes <- names(sort(abs(loadings[,2]), decreasing=TRUE)[1:num_top_genes])

gois <- top_PC1_genes
# gois <- top_PC2_genes

# ---------- Print KEGG descriptions
gois_clean <- sub("\\.[0-9]+$", "", gois)
gene_descriptions <- sapply(gois_clean, function(gene) {
  entry <- tryCatch(
    keggGet(paste0("pfa:", gene)),
    error = function(e) return(NULL)
  )
  if (!is.null(entry) && length(entry) > 0) {
    return(entry[[1]]$NAME)
  } else {
    return(NA)
  }
})

# Create a data frame for easier viewing
desc_df <- data.frame(
  GeneID = gois_clean,
  Description = gene_descriptions
)
print(desc_df)
clean_descs <- sub("(RefSeq) ", "", desc_df[,2], fixed = TRUE)
clean_descs[grepl("unknown function", clean_descs, ignore.case = TRUE)] <- "unknown function"
clean_descs[is.na(clean_descs)] <- "unknown function"

# ---------- PLOT HEATMAP (gene methylation %) of TOP LOADINGS FOR PC1 
annotation_col <- data.frame(
  Group = group,
  Time  = time
)
rownames(annotation_col) <- row_names
annotation_col$Group <- factor(annotation_col$Group)
annotation_col$Time  <- factor(annotation_col$Time)

# Optional colors for annotation
annotation_colors <- list(
  Group = c(control="#1bb6bb", knocksideways="#f46b64"),
  Time  = c("28"="gray", "32"="orange", "36"="green")
)

heatmap_plot <- pheatmap(data[gois,],
         scale="row",               # z-score by gene
         # annotation_col=annotation_col,
         # annotation_colors=annotation_colors,
         show_rownames=TRUE,
         show_colnames=TRUE,
         # main=paste0("PCA1 top ", num_top_genes, " contributing genes (", method, ")"),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         labels_row = clean_descs)

# export 900x1000
ggsave(filename=paste0(name_no_ext, "_heatmap.", ext), plot = heatmap_plot, device=ext, dpi=300, width=9, height=10, units="in")


# ---------- KEGG analysis
# KEGG organism code for Plasmodium falciparum is 'pfa'
# Use gene IDs as they are (PF3D7_XXXXXXX)
# kegg_result <- enrichKEGG(
#   gene         = gois,
#   organism     = 'pfa',
#   keyType      = "kegg", # gene IDs should match KEGG format
#   pvalueCutoff = 0.1
# )

# head(as.data.frame(kegg_result))
# dotplot(kegg_result)

# ---------- GO analysis
gaf <- read.delim(gaf_file, header = FALSE, comment.char = "!", stringsAsFactors = FALSE, sep="\t")
colnames(gaf)[c(2,5,10)] <- c("GeneID", "GO_ID", "Description")
term2gene <- gaf[, c("GO_ID", "GeneID")]
go_ids <- unique(term2gene$GO_ID)
go_names <- Term(go_ids)
term2name <- data.frame(GO_ID = go_ids, GO_Name = go_names)

# go_enrich <- enricher(
#   gene = gois,
#   TERM2GENE = term2gene
# )
# 
# go_ids <- unique(go_enrich@result$ID)
# go_names <- Term(go_ids) # GO.db Term function
# go_term_map <- data.frame(GO_ID = go_ids, GO_Name = go_names)
# 
# go_enrich_df <- as.data.frame(go_enrich@result)
# go_enrich_df <- merge(go_enrich_df, go_term_map, by.x="ID", by.y="GO_ID", all.x=TRUE)
# go_enrich@result$Description <- go_enrich_df$GO_Name
# 
# head(as.data.frame(go_enrich))
# dotplot(go_enrich)

# rank based list
go_gsea <- GSEA(
  geneList   = all_pc1_genes,
  TERM2GENE  = term2gene,
  TERM2NAME  = term2name,
  minGSSize  = 5,
  # maxGSSize = 500,
  pvalueCutoff = 0.01,
  verbose = FALSE
)

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


# ridgeplot(go_gsea)

# go_gsea_df <- as.data.frame(go_gsea)

# cnetplot(go_gsea)

# PC1 transcript count 


# --------- PC1 vs PC2 scatter

library(ggplot2)
library(ggrepel)

pcx <- "PC1"; pcy <- "PC2"
loadings <- as.data.frame(pca$rotation)   # rows = genes, cols = PCs
loadings$gene <- rownames(loadings)

# compute contributions
loadings$contr_PC1 <- (loadings[[pcx]]^2) * 100
loadings$contr_PC2 <- (loadings[[pcy]]^2) * 100
loadings$contr_combined <- loadings$contr_PC1 + loadings$contr_PC2    # percent contribution to PC1+PC2

# choose top genes to label by combined contribution
topLabel <- 20
top_genes <- head(loadings[order(loadings$contr_combined, decreasing = TRUE), "gene"], topLabel)

# top_genes <- "PF3D7_1439600"

# gois_clean <- sub("\\.[0-9]+$", "", top_genes)
# top_genes <- gois_clean
gene_descriptions <- sapply(top_genes, function(gene) {
  entry <- tryCatch(
    keggGet(paste0("pfa:", gene)),
    error = function(e) return(NULL)
  )
  if (!is.null(entry) && length(entry) > 0) {
    return(entry[[1]]$NAME)
  } else {
    return(NA)
  }
})
clean_descs <- sub("(RefSeq) ", "", gene_descriptions, fixed = TRUE)
clean_descs[grepl("unknown function", clean_descs, ignore.case = TRUE)] <- "unknown function"
clean_descs[is.na(clean_descs)] <- "unknown function"
# clean_descs <- stringr::str_wrap(loadings$description, width = 30)

loadings$description <- clean_descs[loadings$gene]


pca_scatter_plot <- ggplot(loadings, aes_string(x = pcx, y = pcy)) +
  geom_hline(yintercept = 0, color = "pink") + geom_vline(xintercept = 0, color = "pink") +
  geom_point(aes(size = contr_combined, color = contr_combined), alpha = 1) +
  scale_size_continuous(name = "Contribution (%)") +
  scale_color_viridis_c(name = "Contribution (%)") +
  geom_text_repel(data = subset(loadings, gene %in% top_genes),
                  aes(label = description), size = 3, max.overlaps = 50) +
  labs(x = paste0(pcx, " (", round(100 * pca$sdev[1]^2 / sum(pca$sdev^2),1), "%)"),
       y = paste0(pcy, " (", round(100 * pca$sdev[2]^2 / sum(pca$sdev^2),1), "%)")) +
       # title = "PCA variable loadings: PC1 vs PC2") +
  theme(legend.position = "none") +
  theme_classic(base_size=20)

ggsave(filename=paste0(name_no_ext, "_pca_scatter.", ext), plot = pca_scatter_plot, device=ext, dpi=300, width=14.8, height=9.6, units="in")

