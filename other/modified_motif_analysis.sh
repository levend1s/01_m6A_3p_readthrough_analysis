#!/bin/bash

set -e

source ./env.sh
# source ${RQC_DIR}/env/bin/activate

OUTDIR=$(pwd)/modified_motif_analysis_output
mkdir -p ${OUTDIR}

# modify plasmodium FASTA based off bedmethyl


READ_DEPTH=15
M6A_RATIO=0.25

# # ---------- NANOPORE BEDMETHYL PROCESSING ---------- #
# FILTER BEDMETHYLS FOR 15X COV AND 25% MOD PROP
# echo "Processing 28hpi samples ..."
# awk -v READ_DEPTH=${READ_DEPTH} -v M6A_RATIO=${M6A_RATIO} '{if ($10 >= READ_DEPTH && ($12/($10+$17)*100) >= M6A_RATIO*100) print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' ${BEDMETHYL_DIR}/28C1_to_pfal_0.95.0p.a.bed > ${OUTDIR}/nanopore_28.temp.bed
# awk -v READ_DEPTH=${READ_DEPTH} -v M6A_RATIO=${M6A_RATIO} '{if ($10 >= READ_DEPTH && ($12/($10+$17)*100) >= M6A_RATIO*100) print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' ${BEDMETHYL_DIR}/28C2_to_pfal_0.95.0p.a.bed >> ${OUTDIR}/nanopore_28.temp.bed
# sort -k1,1 -k2,2n ${OUTDIR}/nanopore_28.temp.bed | uniq > ${OUTDIR}/nanopore_28.bed

# echo "Processing 32hpi samples ..."
# awk -v READ_DEPTH=${READ_DEPTH} -v M6A_RATIO=${M6A_RATIO} '{if ($10 >= READ_DEPTH && ($12/($10+$17)*100) >= M6A_RATIO*100) print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' ${BEDMETHYL_DIR}/32C1_to_pfal_0.95.0p.a.bed > ${OUTDIR}/nanopore_32.temp.bed
# awk -v READ_DEPTH=${READ_DEPTH} -v M6A_RATIO=${M6A_RATIO} '{if ($10 >= READ_DEPTH && ($12/($10+$17)*100) >= M6A_RATIO*100) print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' ${BEDMETHYL_DIR}/32C2_to_pfal_0.95.0p.a.bed >> ${OUTDIR}/nanopore_32.temp.bed
# sort -k1,1 -k2,2n ${OUTDIR}/nanopore_32.temp.bed | uniq > ${OUTDIR}/nanopore_32.bed

# echo "Processing 36hpi samples ..."
# awk -v READ_DEPTH=${READ_DEPTH} -v M6A_RATIO=${M6A_RATIO} '{if ($10 >= READ_DEPTH && ($12/($10+$17)*100) >= M6A_RATIO*100) print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' ${BEDMETHYL_DIR}/36C1_to_pfal_0.95.0p.a.bed > ${OUTDIR}/nanopore_36.temp.bed
# awk -v READ_DEPTH=${READ_DEPTH} -v M6A_RATIO=${M6A_RATIO} '{if ($10 >= READ_DEPTH && ($12/($10+$17)*100) >= M6A_RATIO*100) print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' ${BEDMETHYL_DIR}/36C2_to_pfal_0.95.0p.a.bed >> ${OUTDIR}/nanopore_36.temp.bed
# sort -k1,1 -k2,2n ${OUTDIR}/nanopore_36.temp.bed | uniq > ${OUTDIR}/nanopore_36.bed

# # # FIND UNION
# echo "Finding nanopore union ..."
# cat ${OUTDIR}/nanopore_28.bed > ${OUTDIR}/nanopore_28u32u36.temp.bed
# cat ${OUTDIR}/nanopore_32.bed >> ${OUTDIR}/nanopore_28u32u36.temp.bed
# cat ${OUTDIR}/nanopore_36.bed >> ${OUTDIR}/nanopore_28u32u36.temp.bed
# sort -k1,1 -k2,2n ${OUTDIR}/nanopore_28u32u36.temp.bed | uniq > ${OUTDIR}/nanopore_28u32u36.bed

# ---------- FASTA EDITING ---------- #
# m is negative strand M is forward strand

# echo "Modifying genome FASTA ..."
# awk 'BEGIN{OFS="\t"}
# {
#     pos = $2 + 1         # BED 0-based → VCF 1-based
#     base = ($6=="-") ? "K" : "M"
#     ref = ($6=="-") ? "T" : "A"
#     print $1, pos, ".", ref, base, ".", "PASS", "."
# }' ${OUTDIR}/nanopore_28u32u36.bed > ${OUTDIR}/nanopore_28u32u36.vcf

# printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" \
#   | cat - ${OUTDIR}/nanopore_28u32u36.vcf > ${OUTDIR}/nanopore_28u32u36_with_header.vcf

# bgzip -c ${OUTDIR}/nanopore_28u32u36_with_header.vcf > ${OUTDIR}/nanopore_28u32u36_with_header.vcf.gz
# tabix -p vcf ${OUTDIR}/nanopore_28u32u36_with_header.vcf.gz

# bcftools consensus -f ${GENOME_FILE} ${OUTDIR}/nanopore_28u32u36_with_header.vcf.gz > ${OUTDIR}/plasmodium_genome_modified_motif.fasta

# ---------- TANDEM MOTIF FINDING ---------- #

# rm -f ${OUTDIR}/all_tandem_modified_motif_matches_filtered.bed ${OUTDIR}/all_tandem_unmodified_motif_matches_filtered.bed
# touch ${OUTDIR}/all_tandem_modified_motif_matches_filtered.bed ${OUTDIR}/all_tandem_unmodified_motif_matches_filtered.bed

# for N in {1..30}; do
#   LESS_N=$((N-1))

#   MOTIF=$(printf '%*sM%*sM%*sM%*sM' "$N" '' "$N" '' "$N" '' "$N" '' | tr ' ' '.')
#   echo ${MOTIF}
#   python3 ${RQC_PATH} motif_finder -a ${ANNOTATION_FILE} -g ${OUTDIR}/plasmodium_genome_modified_motif.fasta -o ${OUTDIR}/modified_motif_matches.tsv -m "${MOTIF}"

#   # unmodified motif
#   # UNMODIFIED_MOTIF=$(printf '%*sA[C]%*sA[C]%*sA[C]%*sA[C]' "$LESS_N" '' "$LESS_N" '' "$LESS_N" '' "$LESS_N" '' | tr ' ' '.')
#   # echo ${UNMODIFIED_MOTIF}
#   # python3 ${RQC_PATH} motif_finder -a ${ANNOTATION_FILE} -g ${OUTDIR}/plasmodium_genome_modified_motif.fasta -o ${OUTDIR}/unmodified_motif_matches.tsv -m "${UNMODIFIED_MOTIF}"


#   # filter low complexity matches
#   # COMPLEXITY_THRESHOLD=$(((N+1) * 2))
#   # COMPLEXITY_STRING=$(printf '%*s' "${COMPLEXITY_THRESHOLD}" '' | tr ' ' 'M')
#   COMPLEXITY_STRING="MMMM"
#   tail -n +2 ${OUTDIR}/modified_motif_matches.tsv | grep -v -- "$COMPLEXITY_STRING" | awk '{ sub(/\t+$/, ""); print }' | awk -F'\t' -v OFS='\t' -v N="$N" '{ $5 = N; print }' >> ${OUTDIR}/all_tandem_modified_motif_matches_filtered.bed
  
#   # COMPLEXITY_STRING="AAAA"
#   # tail -n +2 ${OUTDIR}/unmodified_motif_matches.tsv | grep -v -- "$COMPLEXITY_STRING" | awk '{ sub(/\t+$/, ""); print }' | awk -F'\t' -v OFS='\t' -v N="$N" '{ $5 = N; print }' >> ${OUTDIR}/all_tandem_unmodified_motif_matches_filtered.bed

# done

# awk '$3 ~ /^(five_prime_UTR|three_prime_UTR|CDS|ncRNA|rRNA|tRNA|snRNA|snoRNA)$/' ${ANNOTATION_FILE} > ${OUTDIR}/plasmodium_simple_feature.gff # get simple features to a single gff

# bedtools intersect -a ${OUTDIR}/all_tandem_modified_motif_matches_filtered.bed -b ${OUTDIR}/plasmodium_simple_feature.gff -wa -wb -s -loj | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9}' > ${OUTDIR}/modified_matches.bed  # skip 

# sort -k1,1 -k2,2n ${OUTDIR}/modified_matches.bed > ${OUTDIR}/modified_matches.sorted.bed
# bedtools merge -s -c 4,5,6,7 -o first,first,first,first -i ${OUTDIR}/modified_matches.sorted.bed > ${OUTDIR}/modified_matches.sorted.merged.bed

# cat ${OUTDIR}/modified_matches.sorted.merged.bed

# ---------- OR: EXPLICIT MOTIF FINDING ---------- #

# TODO: search for poly A motifs


# # ---------- MOTIF ENRICHMENT ANALYSIS ---------- #

# grab a window 5' of each motif match with respect to strand
# enrich each window for motifs
# bedtools flank \
#   -i ${OUTDIR}/modified_matches.sorted.merged.bed \
#   -g ${GENOME_FAI_FILE} \
#   -l 1000 -r 1000 \
#   -s \
# > ${OUTDIR}/modified_matches.upstream.bed

# bedtools flank \
#   -i ${OUTDIR}/nanopore_28u32u36.bed  \
#   -g ${GENOME_FAI_FILE} \
#   -l 1000 -r 1000 \
#   -s \
# > ${OUTDIR}/modified_matches.upstream.bed

# bedtools getfasta \
#   -fi ${GENOME_FILE} \
#   -bed ${OUTDIR}/modified_matches.upstream.bed \
#   -s \
#   -name \
# > ${OUTDIR}/modified_matches.upstream.fa

# printf "# order 0\nA 0.40\nC 0.10\nG 0.10\nT 0.40\n" > ${OUTDIR}/pfalciparum_upstream.bg
# sed 's/T/U/g' ${OUTDIR}/modified_matches.upstream.fa > ${OUTDIR}/modified_matches.upstream_rna.fa

# # MEME
# ${MEME_PATH} ${OUTDIR}/modified_matches.upstream_rna.fa \
#   -rna \
#   -mod zoops \
#   -minw 4 -maxw 8 \
#   -nmotifs 5 \
#   -bfile ${OUTDIR}/pfalciparum_upstream.bg \
#   -oc ${OUTDIR}/meme_upstream_bg

# grep '<motif ' ${OUTDIR}/meme_upstream_bg/meme.xml | sed -n 's/.*name="\([^"]*\)".*/\1/p' > ${OUTDIR}/motif_names.txt

# --- expand IUPAC codes and find motif locations --- #

# ALTERNATIVE: hardcode a motif to search for
echo "AAAAAA" > ${OUTDIR}/motif_names.txt
echo "UUUUUU" >> ${OUTDIR}/motif_names.txt
echo "WWWWWW" >> ${OUTDIR}/motif_names.txt
echo "SSSSSS" >> ${OUTDIR}/motif_names.txt




# # loop through each line
while read -r seq; do
    # expand IUPAC codes into bracketed options
    expanded_seq="$seq"
    expanded_seq="${expanded_seq//U/T}"

    expanded_seq="${expanded_seq//R/[AG]}"
    expanded_seq="${expanded_seq//Y/[CT]}"
    expanded_seq="${expanded_seq//S/[GC]}"
    expanded_seq="${expanded_seq//K/[GT]}"
    expanded_seq="${expanded_seq//M/[AC]}"
    expanded_seq="${expanded_seq//W/[AT]}"
    expanded_seq="${expanded_seq//B/[CGT]}"
    expanded_seq="${expanded_seq//D/[AGT]}"
    expanded_seq="${expanded_seq//H/[ACT]}"
    expanded_seq="${expanded_seq//V/[ACG]}"
    expanded_seq="${expanded_seq//N/[ACGT]}"

    # do other commands here, e.g., print
    echo "Original: $seq"
    echo "Expanded: $expanded_seq"

    # example: pass to some command
    # my_command "$expanded_seq"
    python3 ${RQC_PATH} motif_finder --disable_lookahead -a ${ANNOTATION_FILE} -g ${GENOME_FILE} -o ${OUTDIR}/${seq}.tsv -m "${expanded_seq}"
    tail -n +2 ${OUTDIR}/${seq}.tsv | awk '{ sub(/\t+$/, ""); print }' > ${OUTDIR}/${seq}.bed
    # tail -n +2 ${OUTDIR}/${seq}.tsv | awk '{ sub(/\t+$/, ""); print }' |bedtools intersect -a - -b ${OUTDIR}/pfal_transcripts.bed -wa -wb > ${OUTDIR}/${seq}_mapped.tsv
    # awk 'BEGIN{OFS="\t"} {print $1, $11, $2, $3, $10, $5, $6}' ${OUTDIR}/${seq}_mapped.tsv > ${OUTDIR}/${seq}_mapped.bed

done < ${OUTDIR}/motif_names.txt

MOTIF_LABELS=$(awk -v OUTDIR=${OUTDIR} '{printf "%s %s/%s.tsv ", $0, OUTDIR, $0}' ${OUTDIR}/motif_names.txt | sed 's/ $//')

# # F: nanopore offset from canonical PAS
awk '$6 != 0 {print $1"\t"$2"\t"$6"\t"$6+1"\t"$3"\t"0"\t"$4}' ${TES_ANALYSIS_DIR}/mRNA_tes_analysis_28hpi_compare.tsv > ${OUTDIR}/28hpi_canonical_pas.tsv
awk 'NR==1 {print "contig\tID\tstart\tend\ttype\tscore\tstrand"; next} {print}' ${OUTDIR}/28hpi_canonical_pas.tsv > temp.tsv && mv temp.tsv ${OUTDIR}/28hpi_canonical_pas.tsv

# # calculate offsets from canonical PAS to motif locations
# python3 ${RQC_PATH} calculate_offsets -s -d 500 -r ${OUTDIR}/28hpi_canonical_pas.tsv -o ${OUTDIR}/motiff_offsets.tsv \
#     -i m6a ${GLORI_X_NANOPORE_DIR}/nanopore_28u32u36.bed ${MOTIF_LABELS}
# python3 ${RQC_PATH} plot_relative_distance -l "canonical m6A" -d 500 -i ${OUTDIR}/motiff_offsets.tsv -o ${OUTDIR}/motiff_offsets.png

python3 ${RQC_PATH} calculate_offsets -s -b -d 500 -r ${GLORI_X_NANOPORE_DIR}/nanopore_28u32u36.bed -o ${OUTDIR}/motiff_offsets.tsv \
    -i ${MOTIF_LABELS}
python3 ${RQC_PATH} plot_relative_distance -l "canonical m6A" -d 500 -i ${OUTDIR}/motiff_offsets.tsv -o ${OUTDIR}/motiff_offsets.png

# I can update this calculate offsets function to search a narrower window, to try and identify sites which have no U rich motif within 50nt or so