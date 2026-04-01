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


# manual
chimeric_matrix <- "/Users/jlevendis/01_m6A_3p_readthrough_analysis/03_featureCounts_output/subread_featureCounts_results_chimeric.txt"
chimeric_columns <- c(1,2,3,4)
transcript_matrix <- "/Users/jlevendis/01_m6A_3p_readthrough_analysis/03_featureCounts_output/subread_featureCounts_results"
transcript_columns <- c(6,7,8,9)

# chimeric_columns <- c(5,6,7,8)
# chimeric_columns <- c(9,10,11,12)
# ---------- chimeric DE (interaction term between treatment AND chimeric/bulk count, figure 5C)

# BASED OFF 4.?gl8.3 https://bioconductor.posit.co/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

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

columns_chimeric <- columns_chimeric[ !grepl("API|MIT", rownames(columns_chimeric), ignore.case = TRUE), , drop = FALSE ]
columns <- columns[ !grepl("API|MIT", rownames(columns), ignore.case = TRUE), , drop = FALSE ]

columns_chimeric <- x_chimeric[,chimeric_columns]
columns <- x[,transcript_columns]
not_chimeric <- columns - columns_chimeric

merged_df <- cbind(columns_chimeric, not_chimeric)
colnames(merged_df) <- c("C1_CHIMERIC", "C2_CHIMERIC", "K1_CHIMERIC", "K2_CHIMERIC", "C1_NOTCHIMERIC", "C2_NOTCHIMERIC", "K1_NOTCHIMERIC", "K2_NOTCHIMERIC")

#TODO THIS SHOULD JUST BE MRNAS!!!!
HasCoverage <- rowSums(columns >= 10) == 4
# ????
# HasBoth <- rowSums(columns_chimeric) > 0 & rowSums(not_chimeric) > 0
# table(HasCoverage, HasBoth)

y <- DGEList(counts=merged_df)
y <- y[HasCoverage,, keep.lib.sizes=FALSE]
TotalLibSize <- 0.5 * y$samples$lib.size[Methylation=="Me"] +
  + 0.5 * y$samples$lib.size[Methylation=="Un"]

y$samples$lib.size <- rep(TotalLibSize, each=2)


merged_df[]


# percent_chimeric <- (columns_chimeric / columns)
# percent_chimeric <- apply(percent_chimeric, c(1,2), function(x) {
#   if(!is.finite(x)) 0 else x
# })




