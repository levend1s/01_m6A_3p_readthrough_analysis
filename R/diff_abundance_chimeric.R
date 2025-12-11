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
chimeric_matrix <- args[1]
chimeric_columns <- as.numeric(strsplit(args[2], ",")[[1]])
transcript_matrix <- args[3]
transcript_columns <- as.numeric(strsplit(args[4], ",")[[1]])
gaf_file <- args[5]
out_file <- args[6]

# ---------- chimeric DE (interaction term between treatment AND chimeric/bulk count, figure 5C)

# export ANNOTATION=~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff      
# awk '$3=="mRNA" {split ($9,x,/[=;]/); print x[2]}' ${ANNOTATION} > pfal_mRNA_list.tsv
# awk 'NR==FNR { gsub(/\..*$/, "", $1); keep[$1]=1; next } NR==1 { print; next } { if($1 in keep) print }' pfal_mRNA_list.tsv ~/rqc/output/8.3_featureCounts_chimeric.txt > ~/rqc/output/8.3_featureCounts_chimeric_mRNAs.txt 

# count_matrix_chimeric <- "~/rqc/output/8.3_featureCounts_chimeric_mRNAs.txt"
# count_matrix <- "~/rqc/output/8.3_featureCounts"

x_chimeric <- read.delim(chimeric_matrix, header=TRUE, row.names="Geneid")
x <- read.delim(transcript_matrix, header=TRUE, skip=1, row.names="Geneid")

# x <- x[rownames(x) %in% rownames(x_chimeric), ]
# x_chimeric <- x_chimeric[ !grepl("API|MIT", rownames(x_chimeric), ignore.case = TRUE), , drop = FALSE ]
# x <- x[ !grepl("API|MIT", rownames(x), ignore.case = TRUE), , drop = FALSE ]

columns_chimeric <- x_chimeric[,chimeric_columns]
columns <- x[,transcript_columns]

columns_chimeric <- columns_chimeric[rowSums(columns_chimeric) > 0,]
columns <- columns[rownames(columns_chimeric), , drop=FALSE]
not_chimeric <- columns - columns_chimeric

merged_df <- cbind(columns_chimeric, not_chimeric)

group <- factor(rep(c("control", "knocksideways"), each=2, times=2))
chimeric_status <- factor(rep(c("chimeric", "not_chimeric"), each=4))

design <- model.matrix(~ group * chimeric_status)
y <- DGEList(counts=merged_df, group=group)

# filter oritinal library by min count, and get norm factors and lib size from original library size
y_bulk <- DGEList(counts = columns)
y_bulk <- calcNormFactors(y_bulk)

y$samples$norm.factors <- rep(y_bulk$samples$norm.factors, times = 2)
y$samples$lib.size <- rep(y_bulk$samples$lib.size, times=2)

y <- estimateDisp(y, design)

# # Negative binomial GLM (edge v2)
fit <- glmQLFit(y, design, robust = TRUE)
qlf_qlm <- glmQLFTest(fit, coef="groupknocksideways:chimeric_statusnot_chimeric")
qlf <- qlf_qlm

df = as.data.frame(topTags(qlf,n=Inf))

## plot!

point_size <- 2
pval_cutoff <- 0.05

df$negLogFDR <- -log10(df$FDR)
df$logFC <- -df$logFC

df$Significant <- "none"
df$Significant[df$logFC > 0 & df$FDR < pval_cutoff] <- "up"
df$Significant[df$logFC < 0 & df$FDR < pval_cutoff] <- "down"

df_not <- df[df$Significant == "none", ]
df_sig <- df[df$Significant != "none", ]

# ------------ write files ! ------------ #

ext <- file_ext(out_file)
name_no_ext <- file_path_sans_ext(out_file)

# md
p <- ggplot() +
  geom_point(data = df_not, aes(x = logCPM, y = logFC), color = "grey", alpha = 1, size=point_size) +
  geom_point(data = df_sig, aes(x = logCPM, y = logFC, color = Significant), alpha = 1, size=point_size) +
  # geom_hline(yintercept = -log10(0.05), color = "red") +
  scale_color_manual(values = c("down" = "blue", "up" = "red", "none" = "gray")) +
  geom_point(alpha = 1.0) +
  theme_classic(base_size = 20) +
  geom_hline(yintercept = 0, color = "red") +
  labs(x = "logCPM", y = "logFC") +
  theme(legend.position = "none")

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

# volcano
p <- ggplot() +
  geom_point(data = df_not, aes(x = logFC, y = negLogFDR), color = "grey", alpha = 1, size = point_size) +
  geom_point(data = df_sig, aes(x = logFC, y = negLogFDR, color = Significant), alpha = 1, size = point_size) +
  # geom_hline(yintercept = -log10(0.05), color = "red") +
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

gene_list <- sign(df$logFC) * df$F
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


