## 4a. Compare BED files between treatments 

```console
bedtools jaccard -a $PROCESSED/jq1/bam_files/As4.1_JQ1_summit_window_pval_smallpval.bed -b $PROCESSED/jq1/bam_files/As4.1_DMSO_summit_window_pval_smallpval.bed
```

intersection    union   jaccard n_intersections
14767038        29749929        0.496372        81996

```console
bedtools jaccard -a $PROCESSED/A485/bam_files/As4.1_A485_summit_window_pval_smallpval.bed -b $PROCESSED/A485/bam_files/As4.1_DMSO_summit_window_pval_smallpval.bed
```

intersection    union   jaccard n_intersections
12075841        29078804        0.41528 70664

```console
bedtools jaccard -a $PROCESSED/A485/bam_files/As4.1_DMSO_summit_window_pval_smallpval.bed -b $PROCESSED/jq1/bam_files/As4.1_DMSO_summit_window_pval_smallpval.bed
```

intersection    union   jaccard n_intersections
25571683        25571683        1       123592

Compare each treated peak set against one another.

```console
bedtools jaccard -a $PROCESSED/H89/bam_files/As4.1_H89_summit_window_pval_smallpval.bed -b $PROCESSED/jq1/bam_files/As4.1_JQ1_summit_window_pval_smallpval.bed
```

intersection    union   jaccard n_intersections
10626708        26591733        0.399625        61325

```console
bedtools jaccard -a $PROCESSED/H89/bam_files/As4.1_H89_summit_window_pval_smallpval.bed -b $PROCESSED/A485/bam_files/As4.1_A485_summit_window_pval_smallpval.bed
```

intersection    union   jaccard n_intersections
9761082 24095037        0.405108        58811

```console
bedtools jaccard -a $PROCESSED/A485/bam_files/As4.1_A485_summit_window_pval_smallpval.bed -b $PROCESSED/jq1/bam_files/As4.1_JQ1_summit_window_pval_smallpval.bed
```

intersection    union   jaccard n_intersections
10663421        23864825        0.446826        63716

## 4b. Compare just the dynamic peaks between each set.

```console
cd $PROCESSED/ 
mkdir comparisons/ && cd comparisons/

export JQ1_DIR=$PROCESSED/jq1/bam_files/
export H89_DIR=$PROCESSED/H89/bam_files/
export A485_DIR=$PROCESSED/A485/bam_files/

LC_COLLATE=C sort -k1,1 -k2,2n ${JQ1_DIR}/dynamic_peaks.bed > JQ1_dynamic_peaks_sorted.bed
LC_COLLATE=C sort -k1,1 -k2,2n ${H89_DIR}/dynamic_peaks.bed > H89_dynamic_peaks_sorted.bed
LC_COLLATE=C sort -k1,1 -k2,2n ${A485_DIR}/dynamic_peaks.bed > A485_dynamic_peaks_sorted.bed

bedtools jaccard -a H89_dynamic_peaks_sorted.bed -b JQ1_dynamic_peaks_sorted.bed
# jaccard: 0.0239274
bedtools jaccard -a H89_dynamic_peaks_sorted.bed -b A485_dynamic_peaks_sorted.bed
# jaccard: 0.048398
bedtools jaccard -a JQ1_dynamic_peaks_sorted.bed -b A485_dynamic_peaks_sorted.bed
# jaccard: 0.0815062
```

## 4c. Also compare the specific dynamic cluster BED files.

```console
export A485_DIR=$PROCESSED/A485/bam_files/FIMO/tomtom/fimo_composites/
LC_COLLATE=C sort -k1,1 -k2,2n ${A485_DIR}/cluster_bed_cluster1.bed6 > A485_increased_dynamic_peaks.bed
LC_COLLATE=C sort -k1,1 -k2,2n ${A485_DIR}/cluster_bed_cluster2.bed6 > A485_decreased_dynamic_peaks.bed

export H89_DIR=$PROCESSED/H89/bam_files/
LC_COLLATE=C sort -k1,1 -k2,2n ${H89_DIR}/cluster_bed_cluster1.bed6 > H89_decreased_dynamic_peaks.bed
LC_COLLATE=C sort -k1,1 -k2,2n ${H89_DIR}/cluster_bed_cluster2.bed6 > H89_increased_dynamic_peaks.bed

export JQ1_DIR=$PROCESSED/jq1/bam_files/FIMO/tomtom/fimo_composites/
LC_COLLATE=C sort -k1,1 -k2,2n ${JQ1_DIR}/cluster_bed_cluster1.bed6 > JQ1_increased_dynamic_peaks.bed
LC_COLLATE=C sort -k1,1 -k2,2n ${JQ1_DIR}/cluster_bed_cluster2.bed6 > JQ1_decreased_dynamic_peaks.bed

bedtools jaccard -a H89_increased_dynamic_peaks.bed -b A485_increased_dynamic_peaks.bed
bedtools jaccard -a H89_increased_dynamic_peaks.bed -b A485_decreased_dynamic_peaks.bed
bedtools jaccard -a H89_decreased_dynamic_peaks.bed -b A485_increased_dynamic_peaks.bed
bedtools jaccard -a H89_decreased_dynamic_peaks.bed -b A485_decreased_dynamic_peaks.bed

bedtools jaccard -a H89_increased_dynamic_peaks.bed -b JQ1_increased_dynamic_peaks.bed
bedtools jaccard -a H89_increased_dynamic_peaks.bed -b JQ1_decreased_dynamic_peaks.bed
bedtools jaccard -a H89_decreased_dynamic_peaks.bed -b JQ1_increased_dynamic_peaks.bed
bedtools jaccard -a H89_decreased_dynamic_peaks.bed -b JQ1_decreased_dynamic_peaks.bed

bedtools jaccard -a JQ1_increased_dynamic_peaks.bed -b A485_increased_dynamic_peaks.bed
bedtools jaccard -a JQ1_increased_dynamic_peaks.bed -b A485_decreased_dynamic_peaks.bed
bedtools jaccard -a JQ1_decreased_dynamic_peaks.bed -b A485_increased_dynamic_peaks.bed
bedtools jaccard -a JQ1_decreased_dynamic_peaks.bed -b A485_decreased_dynamic_peaks.bed
```

Are the binding sites the same between the increased H89 binding sites and the JQ1 decreased region binding sites?

```console
export H89_DIR=$PROCESSED/H89/bam_files/FIMO/tomtom/fimo_composites/main_figure_beds/
export JQ1_DIR=$PROCESSED/jq1/bam_files/FIMO/tomtom/fimo_composites/main_figure_beds/
export A485_DIR=$PROCESSED/A485/bam_files/FIMO/tomtom/fimo_composites/main_figure_beds/

bedtools jaccard -a ${H89_DIR}/Sp9_2M.bed -b ${JQ1_DIR}/Sp9_2M.bed
# jaccard: 0.853942
bedtools jaccard -a ${H89_DIR}/Sp9_fimo_all.bed -b ${JQ1_DIR}/Sp9_fimo_all.bed
# 0.647501

bedtools jaccard -a ${H89_DIR}/Klf1_2M.bed -b ${JQ1_DIR}/Klf1_2M.bed
# jaccard: 1
bedtools jaccard -a ${H89_DIR}/Klf1_fimo_all.bed -b ${JQ1_DIR}/Klf1_fimo_all.bed
# 0.647501

bedtools jaccard -a ${H89_DIR}/Klf1_fimo_all.bed -b ${A485_DIR}/Klf4_fimo_all.bed
# intersection    union   jaccard n_intersections
# 21135826        31202571        0.677374        112242

bedtools jaccard -a ${JQ1_DIR}/Klf1_fimo_all.bed -b ${A485_DIR}/Klf4_fimo_all.bed
# intersection    union   jaccard n_intersections
# 23609555        30399931        0.776632        121400
```








### 17c. Recalculate bigwigs using the BAM files
 - mm10 effective genome size: 2,730,871,774


 - Look at the genes most affected by H89 treatment in the genome browser
 - Look at the regions most affected by JQ1 and A485 treatment in the genome browser

## 18. Compare A485 scRNA-seq overlapped regions to H89/JQ1

```console
cd $PROCESSED/integration/

export H89_DIR="$PROCESSED/integration/"
export JQ1_DIR="$PROCESSED/jq1/integration_analysis/"
export A485_DIR="$PROCESSED/A485/integration/"

awk 'BEGIN {FS=",";OFS="\t"} NR>1 {print $1, $2, $3, $13, $7, $5}' ${H89_DIR}/H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.csv | LC_COLLATE=C sort -k1,1 -k2,2n - > ${H89_DIR}/H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.bed

cat ${H89_DIR}/H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.bed | tr ' ' '\t' > ${H89_DIR}/H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.bed

awk 'BEGIN {FS=",";OFS="\t"} NR>1 {print $1, $2, $3, $14, $7, $5}' ${H89_DIR}/H89down_P300-H3K27Ac_decATAC_Ren1-LO_olp.csv | LC_COLLATE=C sort -k1,1 -k2,2n - > ${H89_DIR}/H89down_P300-H3K27Ac_decATAC_Ren1-LO_olp.bed

awk 'BEGIN {FS=",";OFS="\t"} NR>1 {print $2, $3, $4, $7"-"$20, $15, $6}' ${A485_DIR}/A485_increased-accessibilty-regions_linked_scRNAseq.csv | LC_COLLATE=C sort -k1,1 -k2,2n - > ${A485_DIR}/A485_increased-accessibilty-regions_linked_scRNAseq.bed

awk 'BEGIN {FS=",";OFS="\t"} NR>1 {print $2, $3, $4, $7"-"$20, $15, $6}' ${A485_DIR}/A485_decreased-accessibilty-regions_linked_scRNAseq.csv | LC_COLLATE=C sort -k1,1 -k2,2n - > ${A485_DIR}/A485_decreased-accessibilty-regions_linked_scRNAseq.bed

awk 'BEGIN {FS=",";OFS="\t"} NR>1 {print $2, $3, $4, $7"-"$20, $15, $6}' ${JQ1_DIR}/JQ1_increased-accessibilty-regions_linked_scRNAseq.csv | LC_COLLATE=C sort -k1,1 -k2,2n - > ${JQ1_DIR}/JQ1_increased-accessibilty-regions_linked_scRNAseq.bed

awk 'BEGIN {FS=",";OFS="\t"} NR>1 {print $2, $3, $4, $7"-"$20, $15, $6}' ${JQ1_DIR}/JQ1_decreased-accessibilty-regions_linked_scRNAseq.csv | LC_COLLATE=C sort -k1,1 -k2,2n - > ${JQ1_DIR}/JQ1_decreased-accessibilty-regions_linked_scRNAseq.bed

bedtools summary -i ${H89_DIR}/H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.bed -g /project/gomezlab/rivanna_config/alias/mm10/fasta/default/mm10.chrom.sizes 

bedtools summary -i ${A485_DIR}/A485_increased-accessibilty-regions_linked_scRNAseq.bed -g /project/gomezlab/rivanna_config/alias/mm10/fasta/default/mm10.chrom.sizes 

bedtools intersect -a ${H89_DIR}/H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.bed -b ${A485_DIR}/A485_increased-accessibilty-regions_linked_scRNAseq.bed -loj -filenames

bedtools intersect -a ${H89_DIR}/H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.bed -b ${A485_DIR}/A485_increased-accessibilty-regions_linked_scRNAseq.bed -c
# Asap1: 2 hits
# Tm4sf1: 2 hits
# Tbc1d2: 1 hit

bedtools intersect -a ${H89_DIR}/H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.bed -b ${JQ1_DIR}/JQ1_increased-accessibilty-regions_linked_scRNAseq.bed -c
# Tbc1d2: 1 hit

bedtools intersect -a ${A485_DIR}/A485_increased-accessibilty-regions_linked_scRNAseq.bed -b ${JQ1_DIR}/JQ1_increased-accessibilty-regions_linked_scRNAseq.bed -c | awk -F'\t' '{sum+=$7;}END{print sum;}' 
# 159
wc -l ${A485_DIR}/A485_increased-accessibilty-regions_linked_scRNAseq.bed
# 1839

bedtools jaccard -a ${H89_DIR}/H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.bed -b ${A485_DIR}/A485_increased-accessibilty-regions_linked_scRNAseq.bed
intersection    union   jaccard n_intersections
1010    1069937 0.000943981     5

bedtools jaccard -a ${H89_DIR}/H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.bed -b ${JQ1_DIR}/JQ1_increased-accessibilty-regions_linked_scRNAseq.bed
intersection    union   jaccard n_intersections
202     983870  0.000205312     1

bedtools jaccard -a ${A485_DIR}/A485_increased-accessibilty-regions_linked_scRNAseq.bed -b ${JQ1_DIR}/JQ1_increased-accessibilty-regions_linked_scRNAseq.bed
intersection    union   jaccard n_intersections
30320   650167  0.0466342       159
```

## 19. Compare TF BED files

```console
export H89_DIR=$PROCESSED/H89/bam_files/FIMO/tomtom/fimo_composites/main_figure_beds/
export JQ1_DIR=$PROCESSED/jq1/bam_files/FIMO/tomtom/fimo_composites/main_figure_beds/
export A485_DIR=$PROCESSED/A485/bam_files/FIMO/tomtom/fimo_composites/main_figure_beds/
```

### 19a. SP

```console
bedtools jaccard -a ${H89_DIR}/Sp9_fimo_all.bed -b ${JQ1_DIR}/Sp9_fimo_all.bed > SP_jaccard.txt
bedtools jaccard -a ${H89_DIR}/Sp9_fimo_all.bed -b ${A485_DIR}/Sp2_fimo_all_2023-10-11.bed >> SP_jaccard.txt
bedtools jaccard -a ${JQ1_DIR}/Sp9_fimo_all.bed -b ${A485_DIR}/Sp2_fimo_all_2023-10-11.bed >> SP_jaccard.txt

bedtools jaccard -a ${H89_DIR}/Sp9_2M.bed -b ${JQ1_DIR}/Sp9_2M.bed > SP_2M_jaccard.txt
bedtools jaccard -a ${H89_DIR}/Sp9_2M.bed -b ${A485_DIR}/Sp2_2M_2023-10-11.bed >> SP_2M_jaccard.txt
bedtools jaccard -a ${JQ1_DIR}/Sp9_2M.bed -b ${A485_DIR}/Sp2_2M_2023-10-11.bed >> SP_2M_jaccard.txt
```

### 19b. KLF

```console
bedtools jaccard -a ${H89_DIR}/Klf1_fimo_all.bed -b ${JQ1_DIR}/Klf1_fimo_all.bed > KLF_jaccard.txt
bedtools jaccard -a ${H89_DIR}/Klf1_fimo_all.bed -b ${A485_DIR}/Klf4_fimo_all_2023-10-11.bed >> KLF_jaccard.txt
bedtools jaccard -a ${JQ1_DIR}/Klf1_fimo_all.bed -b ${A485_DIR}/Klf4_fimo_all_2023-10-11.bed >> KLF_jaccard.txt

bedtools jaccard -a ${H89_DIR}/Klf1_2M.bed -b ${JQ1_DIR}/Klf1_2M.bed > KLF_2M_jaccard.txt
bedtools jaccard -a ${H89_DIR}/Klf1_2M.bed -b ${A485_DIR}/Klf4_2M_2023-10-11.bed >> KLF_2M_jaccard.txt
bedtools jaccard -a ${JQ1_DIR}/Klf1_2M.bed -b ${A485_DIR}/Klf4_2M_2023-10-11.bed >> KLF_2M_jaccard.txt
```

### 19c. AP-1

```console
bedtools jaccard -a ${H89_DIR}/Fra2_fimo_all.bed -b ${JQ1_DIR}/Bach2_fimo_all.bed > AP1_jaccard.txt
bedtools jaccard -a ${H89_DIR}/Fra2_fimo_all.bed -b ${A485_DIR}/Fra2_fimo_all_2023-10-11.bed >> AP1_jaccard.txt
bedtools jaccard -a ${JQ1_DIR}/Bach2_fimo_all.bed -b ${A485_DIR}/Fra2_fimo_all_2023-10-11.bed >> AP1_jaccard.txt

bedtools jaccard -a ${H89_DIR}/Fra2_2M.bed -b ${JQ1_DIR}/Bach2_2M.bed > AP1_2M_jaccard.txt
bedtools jaccard -a ${H89_DIR}/Fra2_2M.bed -b ${A485_DIR}/Fra2_2M_2023-10-11.bed >> AP1_2M_jaccard.txt
bedtools jaccard -a ${JQ1_DIR}/Bach2_2M.bed -b ${A485_DIR}/Fra2_2M_2023-10-11.bed >> AP1_2M_jaccard.txt
```

### 19d. ROR

```console
bedtools jaccard -a ${H89_DIR}/Rorb_fimo_all.bed -b ${JQ1_DIR}/Rorb_fimo_all.bed > ROR_jaccard.txt

bedtools jaccard -a ${H89_DIR}/Rorb_2M.bed -b ${JQ1_DIR}/Rorb_2M.bed > ROR_2M_jaccard.txt
```

### 19e. TEAD

```console
bedtools jaccard -a ${JQ1_DIR}/TEAD4_fimo_all.bed -b ${A485_DIR}/TEAD4_fimo_all_2023-10-11.bed >> TEAD_jaccard.txt

bedtools jaccard -a ${JQ1_DIR}/TEAD4_2M.bed -b ${A485_DIR}/TEAD4_2M_2023-10-11.bed > TEAD_2M_jaccard.txt
```

### 19f. RUNX

```console
bedtools jaccard -a ${JQ1_DIR}/RUNX_fimo_all.bed -b ${A485_DIR}/RUNX_fimo_all_2023-10-11.bed >> RUNX_jaccard.txt

bedtools jaccard -a ${JQ1_DIR}/RUNX_2M.bed -b ${A485_DIR}/RUNX_2M_2023-10-11.bed > RUNX_2M_jaccard.txt
```


## 14. Add files of interest to visualize on the UCSC Genome Browser

```console
export PROJECT_WWW="$WWW/genome_browser/trackHub/mm10/"
export WWW_DIR="$PROCESSED/A485/integration/genome_browser/trackHub/mm10/"
cp ${PROJECT_WWW}/As4.1_DMSO_merged.bw ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_H89_merged.bw ${WWW_DIR}
cp ${PROJECT_WWW}/ctrl.bb ${WWW_DIR}
cp ${PROJECT_WWW}/trt.bb ${WWW_DIR}
cp ${PROJECT_WWW}/ATAC_H89_v_DMSO_consensus_peak.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO-v-H89_all_dynamic_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO-v-H89_nondynamic_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO-v-H89_decreased_dynamic_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO-v-H89_increased_dynamic_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO-v-H89_H3K27Ac_consensus_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO-v-H89_H3K27Ac_DMSO_enriched_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO-v-H89_H3K27Ac_H89_enriched_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO-v-H89_P300_consensus_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO-v-H89_P300_DMSO_enriched_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO-v-H89_P300_H89_enriched_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_H89_scRNAseq_RPKM.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/As4.1_DMSO_scRNAseq_RPKM.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H3K27Ac_DMSO_true_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/H3K27Ac_H89_true_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/H3K27Ac_H89_decATAC_2023-05-22.bb ${WWW_DIR}
cp ${PROJECT_WWW}/H3K27Ac_H89_incATAC_2023-05-22.bb ${WWW_DIR}
```
