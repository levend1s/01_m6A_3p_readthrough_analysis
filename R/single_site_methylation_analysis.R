library(edgeR)
library(ggplot2)
library(pheatmap)
library(tools)


args <- commandArgs(trailingOnly = TRUE)

input_file_28C1 <- args[1]
input_file_28C2 <- args[2]
input_file_28K1 <- args[3]
input_file_28K2 <- args[4]

out_file <- args[5]

ext <- file_ext(out_file)
name_no_ext <- file_path_sans_ext(out_file)

# !!! FILTER LOW M6A COUNT / RATIO
M6A_A_THRESHOLD <- 0.00
M6A_COUNT_THRESHOLD <- 20
READ_DEPTH_THRESHOLD <- 10
FDR_CUTOFF <- 0.2
FDR_THRESHOLD_DE <- 0.05
point_size <- 2

cat("Starting script\n")


# score == valid_cov
# fraction_mod = mod / valid_cov
# valid = mod + other_mod + canonical
# my_valid = mod + other_mod + canonical + diff
# my_fraction_mod = mod / my_valid
bedmethyl_header <- c("contig", "start", "end", "code", "score", "strand", "start_pos", "end_pos", "color", "valid_cov", "fraction_mod", "mod", "canonical", "other_mod", "delete", "fail", "diff", "nocall")
feature_map_header <- c("five_prime", "three_prime", "cds", "ncRNA", "rRNA", "tRNA", "snRNA", "snoRNA")

bedmethyl_header <- c(bedmethyl_header, feature_map_header)

df_28C1 <- read.delim(input_file_28C1, header=FALSE, col.names = bedmethyl_header)
df_28C2 <- read.delim(input_file_28C2, header=FALSE, col.names = bedmethyl_header)
df_28K1 <- read.delim(input_file_28K1, header=FALSE, col.names = bedmethyl_header)
df_28K2 <- read.delim(input_file_28K2, header=FALSE, col.names = bedmethyl_header)

cat("Read in bedmethyl files complete!\n")

df_28C1$U_28C1 <- df_28C1$canonical + df_28C1$other_mod + df_28C1$diff# + df_28C1$fail
df_28C2$U_28C2 <- df_28C2$canonical + df_28C2$other_mod + df_28C2$diff# + df_28C2$fail
df_28K1$U_28K1 <- df_28K1$canonical + df_28K1$other_mod + df_28K1$diff# + df_28K1$fail
df_28K2$U_28K2 <- df_28K2$canonical + df_28K2$other_mod + df_28K2$diff# + df_28K2$fail

df_28C1$M_28C1 <- df_28C1$mod
df_28C2$M_28C2 <- df_28C2$mod
df_28K1$M_28K1 <- df_28K1$mod
df_28K2$M_28K2 <- df_28K2$mod

df_28C1_simple <- df_28C1
df_28C2_simple <- df_28C2
df_28K1_simple <- df_28K1
df_28K2_simple <- df_28K2

df_28C1_simple <- df_28C1_simple[,c("contig", "start", "end", "strand", "M_28C1", "U_28C1")]
df_28C2_simple <- df_28C2_simple[,c("contig", "start", "end", "strand", "M_28C2", "U_28C2")]
df_28K1_simple <- df_28K1_simple[,c("contig", "start", "end", "strand", "M_28K1", "U_28K1")]
df_28K2_simple <- df_28K2_simple[,c("contig", "start", "end", "strand", "M_28K2", "U_28K2")]

df_28C1_simple <- df_28C1_simple[(df_28C1_simple$M_28C1 + df_28C1_simple$U_28C1) >= READ_DEPTH_THRESHOLD,]
df_28C2_simple <- df_28C2_simple[(df_28C2_simple$M_28C2 + df_28C2_simple$U_28C2) >= READ_DEPTH_THRESHOLD,]
df_28K1_simple <- df_28K1_simple[(df_28K1_simple$M_28K1 + df_28K1_simple$U_28K1) >= READ_DEPTH_THRESHOLD,]
df_28K2_simple <- df_28K2_simple[(df_28K2_simple$M_28K2 + df_28K2_simple$U_28K2) >= READ_DEPTH_THRESHOLD,]

# !! merge bedmethyls by site
merged_df <- Reduce(function(x, y) merge(x, y, by = c("contig", "start", "end", "strand")), list(df_28C1_simple, df_28C2_simple, df_28K1_simple, df_28K2_simple))
rownames(merged_df) <- paste(merged_df$contig, merged_df$start, merged_df$end, merged_df$strand, sep="_")
merged_df$contig <- NULL
merged_df$start <- NULL
merged_df$end <- NULL
merged_df$strand <- NULL

merged_df <- merged_df[
  (merged_df$M_28C1 >= M6A_COUNT_THRESHOLD &
     merged_df$M_28C2 >= M6A_COUNT_THRESHOLD)
  &
  (  (merged_df$M_28C1 / (merged_df$M_28C1 + merged_df$U_28C1)) >= M6A_A_THRESHOLD &
       (merged_df$M_28C2 / (merged_df$M_28C2 + merged_df$U_28C2)) >= M6A_A_THRESHOLD
  )
  ,
]

cat("Filtering complete!\n")

# !!! site level DE analysis
# sites which are differentially expressed will be excluded from DM analysis
# Note that since sites are from transcripts, expression of this site is complex (related to transcript expression changes too)
merged_df_tot <- data.frame(
  merged_df$M_28C1 + merged_df$U_28C1,
  merged_df$M_28C2 + merged_df$U_28C2,
  merged_df$M_28K1 + merged_df$U_28K1,
  merged_df$M_28K2 + merged_df$U_28K2)
rownames(merged_df_tot) <- rownames(merged_df)
colnames(merged_df_tot) <- c("depth_28C1", "depth_28C2", "depth_28K1", "depth_28K2")

de_group <- factor(rep(c("control", "knocksideways"), each=2))
de_design <- model.matrix(~de_group)
y <- DGEList(counts=merged_df_tot, group=de_group)

y <- calcNormFactors(y)
y <- estimateDisp(y, de_design)

# # Negative binomial GLM (edge v2)
fit <- glmFit(y, de_design)
qlf_nbglm <- glmLRT(fit, coef="de_groupknocksideways")
qlf <- qlf_nbglm

de_res <- topTags(qlf, n=Inf)$table

# plotMD(qlf)


# !!! DM ANALYSIS AND FILTERING OF DIFFERENTIALLY EXPRESSED SITES
de_rownames <- rownames(de_res[de_res$FDR < FDR_THRESHOLD_DE,])
df_merged_not_de <- merged_df[!rownames(merged_df) %in% de_rownames, ]
df_for_dm <- df_merged_not_de

df_for_dm <- merged_df

dm_group <- factor(rep(c("control", "knocksideways"), each=4))
meth_status <- factor(rep(c("methylated", "unmethylated"), times=4))

de_design <- model.matrix(~ dm_group * meth_status)
y <- DGEList(counts=df_for_dm, group=dm_group)

y <- calcNormFactors(y)
y <- estimateDisp(y, de_design)

# # Negative binomial GLM (edge v2)
fit <- glmFit(y, de_design)
qlf_nbglm <- glmLRT(fit, coef="dm_groupknocksideways:meth_statusunmethylated")
qlf <- qlf_nbglm

# Get results table
dm_res <- topTags(qlf, n=Inf)$table

# plotMD(qlf)

# !!! merge DM results with position information (3p, 5p, CDS)
df_gene_positions <- df_28C1[df_28C1$mod >= READ_DEPTH_THRESHOLD,c("contig", "start", "end", "strand", feature_map_header)]
df_gene_positions$position <- apply(df_gene_positions[, feature_map_header, drop = FALSE], 1, function(x) {
  paste(feature_map_header[x == 1], collapse = ",")
})
df_gene_positions <- df_gene_positions[, !(names(df_gene_positions) %in% feature_map_header)]

rownames(df_gene_positions) <- paste(df_gene_positions$contig, df_gene_positions$start, df_gene_positions$end, df_gene_positions$strand, sep="_")

df_gene_positions$contig <- NULL
df_gene_positions$start <- NULL
df_gene_positions$end <- NULL
df_gene_positions$strand <- NULL

# --- final merge
merged_res <- merge(dm_res, df_gene_positions, by = "row.names")

rownames(merged_res) <- merged_res$Row.names
merged_res$Row.names <- NULL
# also merge count information so we can get control m6A/A ratio
merged_res <- merge(merged_res, merged_df, by = "row.names")
rownames(merged_res) <- merged_res$Row.names
merged_res$Row.names <- NULL

merged_res$negLogFDR <- -log10(merged_res$FDR)
merged_res$negLogFC <- -(merged_res$logFC)
merged_res$m6A_A_control_average <- ((merged_res$M_28C1 + merged_res$M_28C2) / (merged_res$M_28C1 + merged_res$M_28C2 + merged_res$U_28C1 + merged_res$U_28C2))
# merged_res$m6A_A_control_average <- merged_res$m6A_A_control_average - 1e-6
merged_res$logit_m6A_A_control_average <- log(merged_res$m6A_A_control_average / (1 - merged_res$m6A_A_control_average))
merged_res <- merged_res[is.finite(merged_res$logit_m6A_A_control_average), ]

# "None" means not mapped to a UTR or CDS, this means it might fall outside of an annotation OR it is from non-mRNA transcripts.
merged_res$position[merged_res$position == "cds"] <- "CDS"
merged_res$position[merged_res$position == "three_prime"] <- "UTR"
merged_res$position[merged_res$position == "five_prime"] <- "UTR"
merged_res$position[merged_res$position == ""] <- "Unmapped"
# merged_res$position[merged_res$position == "snRNA"] <- "ncRNA"
# merged_res$position[merged_res$position == "snoRNA"] <- "ncRNA"
# merged_res$position[merged_res$position == "rRNA"] <- "ncRNA"
# merged_res$position[merged_res$position == "tRNA"] <- "ncRNA"

merged_res$Group <- "None"
merged_res$Group[merged_res$logFC > 0 & merged_res$FDR < FDR_CUTOFF] <- "SIG"
merged_res$Group[merged_res$logFC < 0 & merged_res$FDR < FDR_CUTOFF] <- "SIG"

cor_val <- cor(merged_res$logFC, merged_res$logit_m6A_A_control_average, method = "pearson")
cor_val
nrow(merged_res)







# --------------- END PROCESSING + START PLOTTING --------------- #




# --------------- ALL SITES --------------- #

# # Split data
df_not <- merged_res[merged_res$Group == "None", ]
df_sig <- merged_res[merged_res$Group != "None", ]

all_no_fdr <- ggplot() +
  # geom_point(data = merged_res, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position),
             # alpha = 1, size=2) +
  geom_point(data = merged_res, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position), alpha = 1, size = point_size) +
  geom_hline(yintercept = 0, color = "red") +
  scale_color_viridis_d(option = "D") +
  labs(x = expression(m^6 * "A / A (%)"), y = "log2 FC") +
  theme_classic(base_size = 20) +
  # scale_size_continuous(range = c(2, 10)) +
  scale_color_manual(values = c("rRNA" = "#468e8c", "ncRNA" = "#bd4f77", "UTR" = "#f9e955", "None" = "lightgrey", "CDS" = "#3e0751", "Unmapped" = "#7ac56e", "snoRNA" = "#ea9953", "up" = "red", "down" = "blue")) +
  scale_x_continuous(
    breaks = log(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) /
                   (1 - c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99))),
    labels = paste0(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) * 100, "%")
  )

ggsave(filename=paste0(name_no_ext, "_all.", ext), plot = all_no_fdr, device=ext, dpi=300, width=7, height=4.8, units="in")

all_with_fdr <- ggplot() +
  # geom_point(data = merged_res, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position),
             # alpha = 1, size=2) +
  geom_point(data = df_not, aes(x = logit_m6A_A_control_average, y = negLogFC), color = "gray90", alpha = 1, size = point_size) +
  geom_point(data = df_sig, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position), alpha = 1, size = point_size) +
  geom_hline(yintercept = 0, color = "red") +
  scale_color_viridis_d(option = "D") +
  labs(x = expression(m^6 * "A / A (%)"), y = "log2 FC") +
  theme_classic(base_size = 20) +
  # scale_size_continuous(range = c(2, 10)) +
  scale_color_manual(values = c("rRNA" = "#468e8c", "ncRNA" = "#bd4f77", "UTR" = "#f9e955", "None" = "lightgrey", "CDS" = "#3e0751", "Unmapped" = "#7ac56e", "snoRNA" = "#ea9953", "up" = "red", "down" = "blue")) +
  scale_x_continuous(
    breaks = log(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) /
                   (1 - c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99))),
    labels = paste0(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) * 100, "%")
  )

ggsave(filename=paste0(name_no_ext, "_all_fdr", FDR_CUTOFF, ".", ext), plot = all_with_fdr, device=ext, dpi=300, width=7, height=4.8, units="in")


# --------------- NCRNA --------------- #

merged_res_ncrna <- merged_res[merged_res$position != "CDS" & merged_res$position != "UTR" & merged_res$position != "Unmapped",]
# # Split data
df_not <- merged_res_ncrna[merged_res_ncrna$Group == "None", ]
df_sig <- merged_res_ncrna[merged_res_ncrna$Group != "None", ]



all_no_fdr <- ggplot() +
  # geom_point(data = merged_res, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position),
             # alpha = 1, size=2) +
  geom_point(data = merged_res_ncrna, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position), alpha = 1, size = point_size) +
  geom_hline(yintercept = 0, color = "red") +
  scale_color_viridis_d(option = "D") +
  labs(x = expression(m^6 * "A / A (%)"), y = "log2 FC") +
  theme_classic(base_size = 20) +
  # scale_size_continuous(range = c(2, 10)) +
  scale_color_manual(values = c("rRNA" = "#468e8c", "ncRNA" = "#bd4f77", "UTR" = "#f9e955", "None" = "lightgrey", "CDS" = "#3e0751", "Unmapped" = "#7ac56e", "snoRNA" = "#ea9953", "up" = "red", "down" = "blue")) +
  scale_x_continuous(
    breaks = log(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) /
                   (1 - c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99))),
    labels = paste0(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) * 100, "%")
  )

ggsave(filename=paste0(name_no_ext, "_ncrna.", ext), plot = all_no_fdr, device=ext, dpi=300, width=7, height=4.8, units="in")

all_with_fdr <- ggplot() +
  # geom_point(data = merged_res, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position),
             # alpha = 1, size=2) +
  geom_point(data = df_not, aes(x = logit_m6A_A_control_average, y = negLogFC), color = "gray90", alpha = 1, size = point_size) +
  geom_point(data = df_sig, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position), alpha = 1, size = point_size) +
  geom_hline(yintercept = 0, color = "red") +
  scale_color_viridis_d(option = "D") +
  labs(x = expression(m^6 * "A / A (%)"), y = "log2 FC") +
  theme_classic(base_size = 20) +
  # scale_size_continuous(range = c(2, 10)) +
  scale_color_manual(values = c("rRNA" = "#468e8c", "ncRNA" = "#bd4f77", "UTR" = "#f9e955", "None" = "lightgrey", "CDS" = "#3e0751", "Unmapped" = "#7ac56e", "snoRNA" = "#ea9953", "up" = "red", "down" = "blue")) +
  scale_x_continuous(
    breaks = log(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) /
                   (1 - c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99))),
    labels = paste0(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) * 100, "%")
  )

ggsave(filename=paste0(name_no_ext, "_ncrna_fdr.", FDR_CUTOFF, ".", ext), plot = all_with_fdr, device=ext, dpi=300, width=7, height=4.8, units="in")

# --------------- UNMAPPED --------------- #

# # only unmapped
merged_res_unmapped <- merged_res[merged_res$position == "Unmapped",]
# # Split data
df_not <- merged_res_unmapped[merged_res_ncrna$Group == "None", ]
df_sig <- merged_res_unmapped[merged_res_ncrna$Group != "None", ]



all_no_fdr <- ggplot() +
  # geom_point(data = merged_res, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position),
             # alpha = 1, size=2) +
  geom_point(data = merged_res_unmapped, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position), alpha = 1, size = point_size) +
  geom_hline(yintercept = 0, color = "red") +
  scale_color_viridis_d(option = "D") +
  labs(x = expression(m^6 * "A / A (%)"), y = "log2 FC") +
  theme_classic(base_size = 20) +
  # scale_size_continuous(range = c(2, 10)) +
  scale_color_manual(values = c("rRNA" = "#468e8c", "ncRNA" = "#bd4f77", "UTR" = "#f9e955", "None" = "lightgrey", "CDS" = "#3e0751", "Unmapped" = "#7ac56e", "snoRNA" = "#ea9953", "up" = "red", "down" = "blue")) +
  scale_x_continuous(
    breaks = log(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) /
                   (1 - c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99))),
    labels = paste0(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) * 100, "%")
  )

ggsave(filename=paste0(name_no_ext, "_unmapped.", ext), plot = all_no_fdr, device=ext, dpi=300, width=6, height=4.8, units="in")

all_with_fdr <- ggplot() +
  # geom_point(data = merged_res, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position),
             # alpha = 1, size=2) +
  geom_point(data = df_not, aes(x = logit_m6A_A_control_average, y = negLogFC), color = "gray90", alpha = 1, size = point_size) +
  geom_point(data = df_sig, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position), alpha = 1, size = point_size) +
  geom_hline(yintercept = 0, color = "red") +
  scale_color_viridis_d(option = "D") +
  labs(x = expression(m^6 * "A / A (%)"), y = "log2 FC") +
  theme_classic(base_size = 20) +
  # scale_size_continuous(range = c(2, 10)) +
  scale_color_manual(values = c("rRNA" = "#468e8c", "ncRNA" = "#bd4f77", "UTR" = "#f9e955", "None" = "lightgrey", "CDS" = "#3e0751", "Unmapped" = "#7ac56e", "snoRNA" = "#ea9953", "up" = "red", "down" = "blue")) +
  scale_x_continuous(
    breaks = log(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) /
                   (1 - c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99))),
    labels = paste0(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) * 100, "%")
  )

ggsave(filename=paste0(name_no_ext, "_unmapped_fdr.", FDR_CUTOFF, ".", ext), plot = all_with_fdr, device=ext, dpi=300, width=7, height=4.8, units="in")

# --------------- NCRNA --------------- #

# only mRNAs
merged_res_mrna <- merged_res[merged_res$position == "CDS" | merged_res$position == "UTR",]
# # Split data
df_not <- merged_res_mrna[merged_res_mrna$Group == "None", ]
df_sig <- merged_res_mrna[merged_res_mrna$Group != "None", ]



all_no_fdr <- ggplot() +
  # geom_point(data = merged_res, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position),
             # alpha = 1, size=2) +
  geom_point(data = merged_res_mrna, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position), alpha = 1, size = point_size) +
  geom_hline(yintercept = 0, color = "red") +
  scale_color_viridis_d(option = "D") +
  labs(x = expression(m^6 * "A / A (%)"), y = "log2 FC") +
  theme_classic(base_size = 20) +
  # scale_size_continuous(range = c(2, 10)) +
  scale_color_manual(values = c("rRNA" = "#468e8c", "ncRNA" = "#bd4f77", "UTR" = "#f9e955", "None" = "lightgrey", "CDS" = "#3e0751", "Unmapped" = "#7ac56e", "snoRNA" = "#ea9953", "up" = "red", "down" = "blue")) +
  scale_x_continuous(
    breaks = log(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) /
                   (1 - c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99))),
    labels = paste0(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) * 100, "%")
  )

ggsave(filename=paste0(name_no_ext, "_mrna.", ext), plot = all_no_fdr, device=ext, dpi=300, width=7, height=4.8, units="in")

all_with_fdr <- ggplot() +
  # geom_point(data = merged_res, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position),
             # alpha = 1, size=2) +
  geom_point(data = df_not, aes(x = logit_m6A_A_control_average, y = negLogFC), color = "gray90", alpha = 1, size = point_size) +
  geom_point(data = df_sig, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position), alpha = 1, size = point_size) +
  geom_hline(yintercept = 0, color = "red") +
  scale_color_viridis_d(option = "D") +
  labs(x = expression(m^6 * "A / A (%)"), y = "log2 FC") +
  theme_classic(base_size = 20) +
  # scale_size_continuous(range = c(2, 10)) +
  scale_color_manual(values = c("rRNA" = "#468e8c", "ncRNA" = "#bd4f77", "UTR" = "#f9e955", "None" = "lightgrey", "CDS" = "#3e0751", "Unmapped" = "#7ac56e", "snoRNA" = "#ea9953", "up" = "red", "down" = "blue")) +
  scale_x_continuous(
    breaks = log(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) /
                   (1 - c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99))),
    labels = paste0(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) * 100, "%")
  )

ggsave(filename=paste0(name_no_ext, "_mrna_fdr.", FDR_CUTOFF, ".", ext), plot = all_with_fdr, device=ext, dpi=300, width=7, height=4.8, units="in")



# MD PLOT, COLOURED BY TAG
# ggplot() +
#   geom_point(data = merged_res[merged_res$Group == "SIG",], aes(x = logCPM, y = negLogFC, color = position),
#              alpha = 1, size = point_size) +
#   geom_hline(yintercept = 0, color = "red") +
#   scale_color_viridis_d(option = "D") +
#   labs(x = "logCPM", y = "log2 FC") +
#   theme_classic(base_size = 20)

# # VOLCANO PLOT
# ggplot() +
#   geom_point(data = df_not, aes(x = negLogFC, y = negLogFDR),
#              color = "grey", alpha = 1, size = point_size) +
#   geom_point(data = df_sig, aes(x = negLogFC, y = negLogFDR, color = position),
#              alpha = 1, size = point_size) +
#   geom_hline(yintercept = -log10(0.05), color = "red") +
#   scale_color_manual(values = c("rRNA" = "purple", "ncRNA" = "green", "UTR" = "red", "None" = "grey", "CDS" = "pink", "up" = "red", "down" = "blue")) +
#   labs(x = "logFC", y = "-log10 FDR") +
#   theme_classic(base_size = 20)

# PLOT M6A/A in control against FC (MD PLOT BUT WITH controlM6A/A instead of logCPM)
# expect low m6A/A to have increase in methylation, high m6A/A have decrease in methylation
# would show CDS are low m6A/A and 3pUTR are mixture of low and high, but high% sites are predominantly in UTR
# cds_utr_no_ <- ggplot() +
#   # geom_point(data = merged_res, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position),
#              # alpha = 1, size=2) +
#   geom_point(data = df_not, aes(x = logit_m6A_A_control_average, y = negLogFC), color = "gray90", alpha = 1, size = point_size) +
#   geom_point(data = df_sig, aes(x = logit_m6A_A_control_average, y = negLogFC, color = position), alpha = 1, size = point_size) +
#   geom_hline(yintercept = 0, color = "red") +
#   scale_color_viridis_d(option = "D") +
#   labs(x = expression(m^6 * "A / A (%)"), y = "log2 FC") +
#   theme_classic(base_size = 20) +
#   # scale_size_continuous(range = c(2, 10)) +
#   scale_color_manual(values = c("rRNA" = "#468e8c", "ncRNA" = "#bd4f77", "UTR" = "#f9e955", "None" = "lightgrey", "CDS" = "#3e0751", "Unmapped" = "#7ac56e", "snoRNA" = "#ea9953", "up" = "red", "down" = "blue")) +
#   scale_x_continuous(
#     breaks = log(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) /
#                    (1 - c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99))),
#     labels = paste0(c(0, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99) * 100, "%")
#   )

# ggsave(filename=paste0(name_no_ext, "_cds_utr.", ext), plot = p, device=ext, dpi=300, width=6, height=4.8, units="in")
