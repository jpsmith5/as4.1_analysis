### 12.Prepare and merge signal tracks of interest 

First, generate FASTA files of BED sequences

#### 12a. H3K27ac

```console
cd $PROCESSED/nf-core/chipseq/H3K27ac_DB/

export REFGENIE_FASTA="$PROCESSED/nf-core/chipseq/default/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1.fa"
export NF_FASTA="$PROCESSED/nf-core/chipseq/H89/genome/genome.fa"

# FASTA of the H3K27ac DiffBind H89 enriched peaks
bedtools getfasta -fi ${NF_FASTA} -bed H3K27ac_H89_enriched_peaks.bed > H3K27ac_H89_enriched_peaks.fa
# FASTA of the H3K27ac DiffBind H89 enriched peaks
bedtools getfasta -fi ${NF_FASTA} -bed H3K27ac_DMSO_enriched_peaks.bed > H3K27ac_DMSO_enriched_peaks.fa

# FASTA of the H3K27ac H89 true peak set (overlapped with ATAC-seq peaks)
bedtools getfasta -fi ${NF_FASTA} -bed H3K27ac_H89_true_peaks.bed > H3K27ac_H89_true_peaks.fa
# FASTA of the H3K27ac DMSO true peak set (overlapped with ATAC-seq peaks)
bedtools getfasta -fi ${NF_FASTA} -bed H3K27ac_DMSO_true_peaks.bed > H3K27ac_DMSO_true_peaks.fa

# Remove duplicate FASTA entries
seqkit rmdup -s < H3K27ac_H89_enriched_peaks.fa > tmp.fa && mv tmp.fa H3K27ac_H89_enriched_peaks.fa
seqkit rmdup -s < H3K27ac_DMSO_enriched_peaks.fa > tmp.fa && mv tmp.fa H3K27ac_DMSO_enriched_peaks.fa

seqkit rmdup -s < H3K27ac_H89_true_peaks.fa > tmp.fa && mv tmp.fa H3K27ac_H89_true_peaks.fa
seqkit rmdup -s < H3K27ac_DMSO_true_peaks.fa > tmp.fa && mv tmp.fa H3K27ac_DMSO_true_peaks.fa
```

```console
multiBigwigSummary bins -b ${WWW_DIR}/As4.1_H3K27ac_H89_REP1.bigWig ${WWW_DIR}/As4.1_H3K27ac_H89_REP2.bigWig --labels H3K27ac_H89_1 H3K27ac_H89_2 -out H3K27ac_H89_scores_per_bin.npz --outRawCounts H3K27ac_H89_scores_per_bin.tab
plotCorrelation \
-in H3K27ac_H89_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H3K27ac_H89_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix H3K27ac_H89_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/As4.1_H3K27ac_H89_REP1.bigWig --bigwig2 ${WWW_DIR}/As4.1_H3K27ac_H89_REP2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/As4.1_H3K27ac_H89.bigWig --outFileFormat bigwig

multiBigwigSummary bins -b ${WWW_DIR}/As4.1_H3K27ac_DMSO_REP1.bigWig ${WWW_DIR}/As4.1_H3K27ac_DMSO_REP2.bigWig --labels H3K27ac_DMSO_1 H3K27ac_DMSO_2 -out H3K27ac_DMSO_scores_per_bin.npz --outRawCounts H3K27ac_DMSO_scores_per_bin.tab
plotCorrelation \
-in H3K27ac_DMSO_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H3K27ac_DMSO_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix H3K27ac_DMSO_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/As4.1_H3K27ac_DMSO_REP1.bigWig --bigwig2 ${WWW_DIR}/As4.1_H3K27ac_DMSO_REP2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/As4.1_H3K27ac_DMSO.bigWig --outFileFormat bigwig
```

#### 12b. P300
```console
cd $PROCESSED/nf-core/chipseq/P300_DB/

export REFGENIE_FASTA="$PROCESSED/nf-core/chipseq/default/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1.fa"
export NF_FASTA="$PROCESSED/nf-core/chipseq/H89/genome/genome.fa"

# FASTA of the P300 DiffBind H89 enriched peaks
bedtools getfasta -fi ${NF_FASTA} -bed P300_H89_enriched_peaks.bed > P300_H89_enriched_peaks.fa
# FASTA of the P300 DiffBind H89 enriched peaks
bedtools getfasta -fi ${NF_FASTA} -bed P300_DMSO_enriched_peaks.bed > P300_DMSO_enriched_peaks.fa

# FASTA of the P300 H89 true peak set (overlapped with ATAC-seq peaks)
bedtools getfasta -fi ${NF_FASTA} -bed P300_H89_true_peaks.bed > P300_H89_true_peaks.fa
# FASTA of the P300 DMSO true peak set (overlapped with ATAC-seq peaks)
bedtools getfasta -fi ${NF_FASTA} -bed P300_DMSO_true_peaks.bed > P300_DMSO_true_peaks.fa

# Remove duplicate FASTA entries
seqkit rmdup -s < P300_H89_enriched_peaks.fa > tmp.fa && mv tmp.fa P300_H89_enriched_peaks.fa
seqkit rmdup -s < P300_DMSO_enriched_peaks.fa > tmp.fa && mv tmp.fa P300_DMSO_enriched_peaks.fa

seqkit rmdup -s < P300_H89_true_peaks.fa > tmp.fa && mv tmp.fa P300_H89_true_peaks.fa
seqkit rmdup -s < P300_DMSO_true_peaks.fa > tmp.fa && mv tmp.fa P300_DMSO_true_peaks.fa
```

```console
multiBigwigSummary bins -b ${WWW_DIR}/As4.1_P300_H89_REP1.bigWig ${WWW_DIR}/As4.1_P300_H89_REP2.bigWig --labels P300_H89_1 P300_H89_2 -out P300_H89_scores_per_bin.npz --outRawCounts P300_H89_scores_per_bin.tab
plotCorrelation \
-in P300_H89_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o P300_H89_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix P300_H89_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/As4.1_P300_H89_REP1.bigWig --bigwig2 ${WWW_DIR}/As4.1_P300_H89_REP2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/As4.1_P300_H89.bigWig --outFileFormat bigwig

multiBigwigSummary bins -b ${WWW_DIR}/As4.1_P300_DMSO_REP1.bigWig ${WWW_DIR}/As4.1_P300_DMSO_REP2.bigWig --labels P300_DMSO_1 P300_DMSO_2 -out P300_DMSO_scores_per_bin.npz --outRawCounts P300_DMSO_scores_per_bin.tab
plotCorrelation \
-in P300_DMSO_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o P300_DMSO_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix P300_DMSO_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/As4.1_P300_DMSO_REP1.bigWig --bigwig2 ${WWW_DIR}/As4.1_P300_DMSO_REP2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/As4.1_P300_DMSO.bigWig --outFileFormat bigwig
```

#### 12c. Move signal tracks to genome_browser/ location

```console
export H3K27ac_BW="$PROCESSED/chipseq/H89/bwa/mergedLibrary/bigwig/"
cp ${H3K27ac_BW}/As4.1_H3K27ac_DMSO_REP1.bigWig ${WWW_DIR}
cp ${H3K27ac_BW}/As4.1_H3K27ac_DMSO_REP2.bigWig ${WWW_DIR}
cp ${H3K27ac_BW}/As4.1_H3K27ac_H89_REP1.bigWig ${WWW_DIR}
cp ${H3K27ac_BW}/As4.1_H3K27ac_H89_REP2.bigWig ${WWW_DIR}
cp ${H3K27ac_BW}/As4.1_H3K27ac_pooled.bigWig ${WWW_DIR}/As4.1_ChIP_pooled.bigWig

export P300_BW="$PROCESSED/chipseq/P300/bwa/mergedLibrary/bigwig/"
cp ${P300_BW}/As4.1_P300_DMSO_REP1.bigWig ${WWW_DIR}
cp ${P300_BW}/As4.1_P300_DMSO_REP2.bigWig ${WWW_DIR}
cp ${P300_BW}/As4.1_P300_H89_REP1.bigWig ${WWW_DIR}
cp ${P300_BW}/As4.1_P300_H89_REP2.bigWig ${WWW_DIR}
```