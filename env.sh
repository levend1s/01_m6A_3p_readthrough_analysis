#!/bin/bash

# USER SPECIFIED
FEATURECOUNTS_PATH="/Users/jlevendis/Downloads/subread-2.1.1-macOS-arm64/bin/featureCounts"
MODKIT_PATH="/Users/jlevendis/modkit/target/debug/modkit"
RQC_DIR="/Users/jlevendis/rqc"
RQC_PATH="${RQC_DIR}/rqc.py"

PLASMODB_DIR="/Users/jlevendis/Downloads/Pfalciparum3D7/"
ANNOTATION_FILE="${PLASMODB_DIR}/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"
GENOME_FILE="${PLASMODB_DIR}/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"
GENOME_FAI_FILE="${PLASMODB_DIR}/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta.fai"
GAF_FILE="${PLASMODB_DIR}/gaf/PlasmoDB-67_Pfalciparum3D7_GO.gaf"

BAMDIR="/Users/jlevendis/Downloads/sra_upload_bams"

FILTERED_BAMDIR=$(pwd)/01_BAM_filtering_output

# requirements: 28C1_to_pfal_0.95.0p.a.bed, 28C2_to_pfal_0.95.0p.a.bed, 32C1_to_pfal_0.95.0p.a.bed, 32C2_to_pfal_0.95.0p.a.bed, 36C1_to_pfal_0.95.0p.a.bed, 36C2_to_pfal_0.95.0p.a.bed
BEDMETHYL_DIR=$(pwd)/02_generating_bedmethyls_output
FEATURECOUNTS_DIR=$(pwd)/03_featureCounts_output

TES_ANALYSIS_DIR=$(pwd)/04_tes_and_gene_methylation_analysis_output
GLORI_X_NANOPORE_DIR=$(pwd)/05_glori_x_nanopore_overlap_output

GENE_NEIGHBOUR_DIR=$(pwd)/06_gene_neighbours_output
FEATURE_MAPPED_BEDMETHYL_DIR=$(pwd)/07_bedmethyl_feature_mapping_output

# requirements: glori-seq_12hpi.tsv, glori-seq_24hpi.tsv, glori-seq_48hpi.tsv
GLORI_DIR=/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/rqc_output/GLORI-seq

