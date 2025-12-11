# install.packages("eulerr")

library(eulerr)
library(ggplot2)
library(tools)
library(ggplotify)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
glori_x_nanopore_dir <- args[1]
out_file <- args[2]


ext <- file_ext(out_file)
name_no_ext <- file_path_sans_ext(out_file)

# Helper to run shell commands and return numeric
run_int <- function(cmd) {
  out <- system(cmd, intern = TRUE)
  as.numeric(out)
}

# Thanks ChatGPT for turning the 04_glori_x_nanopore_overlap.sh logic into the following R code
# ---------- TIMEPOINT COUNTS ---------- #
cmd_glori_12hpi_count <- sprintf("wc -l < %s/glori_12.bed", glori_x_nanopore_dir)
cmd_glori_24hpi_count <- sprintf("wc -l < %s/glori_24.bed", glori_x_nanopore_dir)
cmd_glori_48hpi_count <- sprintf("wc -l < %s/glori_48.bed", glori_x_nanopore_dir)

cmd_nanopore_28hpi_count <- sprintf("wc -l < %s/nanopore_28.bed", glori_x_nanopore_dir)
cmd_nanopore_32hpi_count <- sprintf("wc -l < %s/nanopore_32.bed", glori_x_nanopore_dir)
cmd_nanopore_36hpi_count <- sprintf("wc -l < %s/nanopore_36.bed", glori_x_nanopore_dir)


# ---------- NANOPORE VENN ---------- #
cmd_np_28_only <- sprintf(
  "bedtools intersect -v -s -a %s/nanopore_28.bed -b %s/nanopore_32.bed | bedtools intersect -v -s -a - -b %s/nanopore_36.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_np_32_only <- sprintf(
  "bedtools intersect -v -s -a %s/nanopore_32.bed -b %s/nanopore_28.bed | bedtools intersect -v -s -a - -b %s/nanopore_36.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_np_36_only <- sprintf(
  "bedtools intersect -v -s -a %s/nanopore_36.bed -b %s/nanopore_32.bed | bedtools intersect -v -s -a - -b %s/nanopore_28.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_np_28_32_only <- sprintf(
  "bedtools intersect -u -s -a %s/nanopore_28.bed -b %s/nanopore_32.bed | bedtools intersect -v -s -a - -b %s/nanopore_36.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_np_28_36_only <- sprintf(
  "bedtools intersect -u -s -a %s/nanopore_28.bed -b %s/nanopore_36.bed | bedtools intersect -v -s -a - -b %s/nanopore_32.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_np_32_36_only <- sprintf(
  "bedtools intersect -u -s -a %s/nanopore_32.bed -b %s/nanopore_36.bed | bedtools intersect -v -s -a - -b %s/nanopore_28.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_np_28_32_36_all <- sprintf(
  "bedtools intersect -u -s -a %s/nanopore_28.bed -b %s/nanopore_32.bed | bedtools intersect -u -s -a - -b %s/nanopore_36.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)


# ---------- GLORI-seq VENN ---------- #
cmd_glori_12_only <- sprintf(
  "bedtools intersect -v -s -a %s/glori_12.bed -b %s/glori_24.bed | bedtools intersect -v -s -a - -b %s/glori_48.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_glori_24_only <- sprintf(
  "bedtools intersect -v -s -a %s/glori_24.bed -b %s/glori_12.bed | bedtools intersect -v -s -a - -b %s/glori_48.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_glori_48_only <- sprintf(
  "bedtools intersect -v -s -a %s/glori_48.bed -b %s/glori_12.bed | bedtools intersect -v -s -a - -b %s/glori_24.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_glori_12_24_only <- sprintf(
  "bedtools intersect -u -s -a %s/glori_12.bed -b %s/glori_24.bed | bedtools intersect -v -s -a - -b %s/glori_48.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_glori_12_48_only <- sprintf(
  "bedtools intersect -u -s -a %s/glori_12.bed -b %s/glori_48.bed | bedtools intersect -v -s -a - -b %s/glori_24.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_glori_24_48_only <- sprintf(
  "bedtools intersect -u -s -a %s/glori_24.bed -b %s/glori_48.bed | bedtools intersect -v -s -a - -b %s/glori_12.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_glori_12_24_48_all <- sprintf(
  "bedtools intersect -u -s -a %s/glori_12.bed -b %s/glori_24.bed | bedtools intersect -u -s -a - -b %s/glori_48.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir, glori_x_nanopore_dir)


# ---------- TOTAL DETECTED ACROSS ALL TIMEPOINTS ---------- #
cmd_np_union_count <- sprintf("wc -l < %s/nanopore_28u32u36.bed", glori_x_nanopore_dir)
cmd_glori_union_count <- sprintf("wc -l < %s/glori_12u24u48.bed", glori_x_nanopore_dir)

cmd_np_only <- sprintf(
  "bedtools intersect -v -s -a %s/nanopore_28u32u36.bed -b %s/glori_12u24u48.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_glori_only <- sprintf(
  "bedtools intersect -v -s -a %s/glori_12u24u48.bed -b %s/nanopore_28u32u36.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_glori_np_union_intersect <- sprintf("wc -l < %s/glori_nanopore_union_intersection.bed", glori_x_nanopore_dir)


# ---------- ONLY SITES PRESENT AT ALL TIMEPOINTS ---------- #
cmd_np_intersection_count <- sprintf("wc -l < %s/nanopore_28i32i36.bed", glori_x_nanopore_dir)
cmd_glori_intersection_count <- sprintf("wc -l < %s/glori_12i24i48.bed", glori_x_nanopore_dir)

cmd_np_intersection_only <- sprintf(
  "bedtools intersect -v -s -a %s/nanopore_28i32i36.bed -b %s/glori_12i24i48.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_glori_intersection_only <- sprintf(
  "bedtools intersect -v -s -a %s/glori_12i24i48.bed -b %s/nanopore_28i32i36.bed | wc -l",
  glori_x_nanopore_dir, glori_x_nanopore_dir)

cmd_glori_np_intersection_intersect <- sprintf(
  "wc -l < %s/glori_nanopore_intersection_intersection.bed",
  glori_x_nanopore_dir)

# ---------- TOTAL DETECTED AT ALL TIMEPOINTS SPLIT BY DRACHS ---------- #
cmd_np_drach_count <- sprintf("wc -l < %s/nanopore_drach_m6A_sites_28u32u36.bed", glori_x_nanopore_dir)
cmd_np_non_drach_count <- sprintf("wc -l < %s/nanopore_non_drach_m6A_sites_28u32u36.bed", glori_x_nanopore_dir)

cmd_glori_drach_count <- sprintf("wc -l < %s/glori_drach_m6A_sites_12u24u48.bed", glori_x_nanopore_dir)
cmd_glori_non_drach_count <- sprintf("wc -l < %s/glori_non_drach_m6A_sites_12u24u48.bed", glori_x_nanopore_dir)

# ------------ # 

sites_intersect_all_timepoints <- c(
    "ONT" = run_int(cmd_np_intersection_only),
    "GLORI" = run_int(cmd_glori_intersection_only),
    "ONT&GLORI" = run_int(cmd_glori_np_intersection_intersect)
)

fit <- euler(sites_intersect_all_timepoints, shape = "ellipse")
p <- as.ggplot(plot(
  fit, 
  fills=c("light blue", "orange"), 
  alpha=1.0, 
  edges=FALSE, 
  # labels = list(cex = 2),
  quantities = list(cex = 1.5),
  labels=FALSE,
  # quantities = TRUE
))

# Save using ggsave
ggsave(
  filename = paste0(name_no_ext, "_intersect_all_timepoints.png"),
  plot = p,
  width = 6,
  height = 6,
  dpi = 300
)

# ---------------------------------- #

sites_total_all_timepoints <- c(
    "ONT" = run_int(cmd_np_only),
    "GLORI" = run_int(cmd_glori_only),
    "ONT&GLORI" = run_int(cmd_glori_np_union_intersect)
)

fit <- euler(sites_total_all_timepoints, shape = "ellipse")
p <- as.ggplot(plot(
  fit, 
  fills=c("light blue", "orange"), 
  alpha=1.0, 
  edges=FALSE, 
  # labels = list(cex = 2),
  quantities = list(cex = 1.5),
  labels=FALSE,
  # quantities = TRUE
))

# Save using ggsave
ggsave(
  filename = paste0(name_no_ext, "_total_all_timepoints.png"),
  plot = p,
  width = 6,
  height = 6,
  dpi = 300
)

# ---------------------------------- #

# Build the vector dynamically
sites_nanopore <- c(
  "28hpi" = run_int(cmd_np_28_only),
  "32hpi" = run_int(cmd_np_32_only),
  "36hpi" = run_int(cmd_np_36_only),
  "28hpi&32hpi" = run_int(cmd_np_28_32_only),
  "28hpi&36hpi" = run_int(cmd_np_28_36_only),
  "32hpi&36hpi" = run_int(cmd_np_32_36_only),
  "28hpi&32hpi&36hpi" = run_int(cmd_np_28_32_36_all)
)

fit <- euler(sites_nanopore, shape = "ellipse")
p <- as.ggplot(plot(
  fit, 
  fills=c("pink", "darkslategray2", "bisque"), 
  alpha=1.0, 
  edges=FALSE, 
  # labels = list(cex = 2),
  quantities = list(cex = 1.5),
  labels=FALSE,
  # quantities = TRUE
))

# Save using ggsave
ggsave(
  filename = paste0(name_no_ext, "_sites_nanopore.png"),
  plot = p,
  width = 6,
  height = 6,
  dpi = 300
)

# ---------------------------------- #


sites_glori <- c(
  "12hpi" = run_int(cmd_glori_12_only),
  "24hpi" = run_int(cmd_glori_24_only),
  "48hpi" = run_int(cmd_glori_48_only),
  "12hpi&24hpi" = run_int(cmd_glori_12_24_only),
  "24hpi&48hpi" = run_int(cmd_glori_24_48_only),
  "12hpi&48hpi" = run_int(cmd_glori_12_48_only),
  "12hpi&24hpi&48hpi" = run_int(cmd_glori_12_24_48_all)
)

fit <- euler(sites_glori, shape = "ellipse")
p <- as.ggplot(plot(
  fit, 
  fills=c("pink", "darkslategray2", "bisque"), 
  alpha=1.0, 
  edges=FALSE, 
  # labels = list(cex = 2),
  quantities = list(cex = 1.5),
  labels=FALSE,
  # quantities = TRUE
))

# Save using ggsave
ggsave(
  filename = paste0(name_no_ext, "_sites_glori.png"),
  plot = p,
  width = 6,
  height = 6,
  dpi = 300
)


# ---------------------------------- #



# nanopore
df <- data.frame(
  category = c("non-DRACH", "DRACH"),
  value = c(run_int(cmd_np_non_drach_count), run_int(cmd_np_drach_count))
)

df <- df %>%
  mutate(
    pct = round(value / sum(value) * 100),
    label = paste0(category, " \n", value, " (", pct, "%)"),
    ypos = cumsum(value) - 0.5*value
  )

# Pie chart
p <- ggplot(df, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(y = ypos, label = label), color = "white", size=6) +
  theme_void(base_size=20) + 
  theme(legend.position = "none")

ggsave(
  filename = paste0(name_no_ext, "_sites_nanopore_non_drach.png"),
  plot = p,
  width = 6,
  height = 6,
  dpi = 300
)

# ---------------------------------- #

# glori
df <- data.frame(
  category = c("non-DRACH", "DRACH"),
  value = c(run_int(cmd_glori_non_drach_count), run_int(cmd_glori_drach_count))
)

df <- df %>%
  mutate(
    pct = round(value / sum(value) * 100),
    label = paste0(category, " \n", value, " (", pct, "%)"),
    ypos = cumsum(value) - 0.5*value
  )

# Pie chart
p <- ggplot(df, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(y = ypos, label = label), color = "white", size=6) +
  theme_void(base_size=20) + 
  theme(legend.position = "none")

ggsave(
  filename = paste0(name_no_ext, "_sites_glori_non_drach.png"),
  plot = p,
  width = 6,
  height = 6,
  dpi = 300
)