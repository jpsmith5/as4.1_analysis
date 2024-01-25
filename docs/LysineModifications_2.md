## 2 Compare signal tracks

### 2a. H3K27Ac ChIP-seq

```console
multiBigwigSummary bins -b ${WWW_DIR}/As4.1_H3K27Ac_H89_REP1.bigWig ${WWW_DIR}/As4.1_H3K27Ac_H89_REP2.bigWig --labels H3K27Ac_H89_1 H3K27Ac_H89_2 -out H3K27Ac_H89_scores_per_bin.npz --outRawCounts H3K27Ac_H89_scores_per_bin.tab
plotCorrelation \
-in H3K27Ac_H89_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H3K27Ac_H89_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix H3K27Ac_H89_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/As4.1_H3K27Ac_H89_REP1.bigWig --bigwig2 ${WWW_DIR}/As4.1_H3K27Ac_H89_REP2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/As4.1_H3K27Ac_H89.bigWig --outFileFormat bigwig

multiBigwigSummary bins -b ${WWW_DIR}/As4.1_H3K27Ac_DMSO_REP1.bigWig ${WWW_DIR}/As4.1_H3K27Ac_DMSO_REP2.bigWig --labels H3K27Ac_DMSO_1 H3K27Ac_DMSO_2 -out H3K27Ac_DMSO_scores_per_bin.npz --outRawCounts H3K27Ac_DMSO_scores_per_bin.tab
plotCorrelation \
-in H3K27Ac_DMSO_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H3K27Ac_DMSO_PearsonCorr_bigwigScores.png   \
--outFileCorMatrix H3K27Ac_DMSO_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/As4.1_H3K27Ac_DMSO_REP1.bigWig --bigwig2 ${WWW_DIR}/As4.1_H3K27Ac_DMSO_REP2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/As4.1_H3K27Ac_DMSO.bigWig --outFileFormat bigwig
```

### 2b. Generate comparative heatmaps of the lysine modifications in untreated cells

 - H3K27Ac
 - H3K4me1
 - H3K4me3
 - H4K16ac
 - H4K5ac
 - H2B5Kac

#### H3K27Ac Cut n Tag

Now, check across all variants of this mark, both the Cut and Tag and the ChIP-seq, then combine
```console
# Cut and Tag H3K27Ac
multiBigwigSummary bins -b ${WWW_DIR}/H3K27ac_R1.bigWig \
 ${WWW_DIR}/H3K27ac_R2.bigWig \
 --labels H3K27Ac_CutNTag_1 H3K27Ac_CutNTag_2 \
 -out H3K27Ac_CutNTag_scores_per_bin.npz \
 --outRawCounts H3K27Ac_CutNTag_scores_per_bin.tab

plotCorrelation \
-in H3K27Ac_CutNTag_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H3K27Ac_CutNTag_scatter_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix H3K27Ac_CutNTag_scatter_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/H3K27ac_R1.bigWig --bigwig2 ${WWW_DIR}/H3K27ac_R2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/H3K27ac_CutNTag.bigWig --outFileFormat bigwig

# All together
multiBigwigSummary bins -b ${WWW_DIR}/As4.1_H3K27Ac_DMSO_REP1.bigWig \
 ${WWW_DIR}/As4.1_H3K27Ac_DMSO_REP2.bigWig \
 ${WWW_DIR}/H3K27ac_R1.bigWig \
 ${WWW_DIR}/H3K27ac_R2.bigWig \
 --labels H3K27Ac_DMSO_1 H3K27Ac_DMSO_2 H3K27Ac_CutNTag_1 H3K27Ac_CutNTag_2 \
 -out H3K27Ac_As4.1_scores_per_bin.npz \
 --outRawCounts H3K27Ac_As4.1_scores_per_bin.tab

plotCorrelation \
-in H3K27Ac_As4.1_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H3K27Ac_As4.1_scatter_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix H3K27Ac_As4.1_scatter_PearsonCorr_bigwigScores.tab

plotCorrelation \
-in H3K27Ac_As4.1_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H3K27Ac_As4.1_scatter_SpearmanCorr_bigwigScores.svg   \
--outFileCorMatrix H3K27Ac_As4.1_scatter_SpearmanCorr_bigwigScores.tab

plotCorrelation \
-in H3K27Ac_As4.1_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o H3K27Ac_As4.1_heatmap_PearsonCorr_bigwigScores.svg \
--outFileCorMatrix H3K27Ac_As4.1_heatmap_PearsonCorr_bigwigScores.tab

plotCorrelation \
-in H3K27Ac_As4.1_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Average Scores Per Bin" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o H3K27Ac_As4.1_heatmap_SpearmanCorr_bigwigScores.svg \
--outFileCorMatrix H3K27Ac_As4.1_heatmap_SpearmanCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/H3K27ac_CutNTag.bigWig --bigwig2 ${WWW_DIR}/As4.1_H3K27Ac_DMSO.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/H3K27Ac_As4.1.bigWig --outFileFormat bigwig
```

#### H3K4me1

```console
multiBigwigSummary bins -b ${WWW_DIR}/H3K4me1_R1.bigWig \
 ${WWW_DIR}/H3K4me1_R2.bigWig \
 --labels H3K4me1_CutNTag_1 H3K4me1_CutNTag_2 \
 -out H3K4me1_scores_per_bin.npz \
 --outRawCounts H3K4me1_scores_per_bin.tab

plotCorrelation \
-in H3K4me1_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H3K4me1_scatter_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix H3K4me1_scatter_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/H3K4me1_R1.bigWig --bigwig2 ${WWW_DIR}/H3K4me1_R2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/H3K4me1.bigWig --outFileFormat bigwig
```

#### H3K4me3

```console
multiBigwigSummary bins -b ${WWW_DIR}/H3K4me3_R1.bigWig \
 ${WWW_DIR}/H3K4me3_R2.bigWig \
 --labels H3K4me3_CutNTag_1 H3K4me3_CutNTag_2 \
 -out H3K4me3_scores_per_bin.npz \
 --outRawCounts H3K4me3_scores_per_bin.tab

plotCorrelation \
-in H3K4me3_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H3K4me3_scatter_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix H3K4me3_scatter_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/H3K4me3_R1.bigWig --bigwig2 ${WWW_DIR}/H3K4me3_R2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/H3K4me3.bigWig --outFileFormat bigwig
```

#### H4K16ac

```console
multiBigwigSummary bins -b ${WWW_DIR}/H4K16ac_R1.bigWig \
 ${WWW_DIR}/H4K16ac_R2.bigWig \
 --labels H4K16ac_CutNTag_1 H4K16ac_CutNTag_2 \
 -out H4K16ac_scores_per_bin.npz \
 --outRawCounts H4K16ac_scores_per_bin.tab

plotCorrelation \
-in H4K16ac_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H4K16ac_scatter_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix H4K16ac_scatter_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/H4K16ac_R1.bigWig --bigwig2 ${WWW_DIR}/H4K16ac_R2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/H4K16ac.bigWig --outFileFormat bigwig
```

#### H4K5ac

```console
multiBigwigSummary bins -b ${WWW_DIR}/H4K5ac_R1.bigWig \
 ${WWW_DIR}/H4K5ac_R2.bigWig \
 --labels H4K5ac_CutNTag_1 H4K5ac_CutNTag_2 \
 -out H4K5ac_scores_per_bin.npz \
 --outRawCounts H4K5ac_scores_per_bin.tab

plotCorrelation \
-in H4K5ac_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H4K5ac_scatter_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix H4K5ac_scatter_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/H4K5ac_R1.bigWig --bigwig2 ${WWW_DIR}/H4K5ac_R2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/H4K5ac.bigWig --outFileFormat bigwig
```

#### H2B5Kac

```console
multiBigwigSummary bins -b ${WWW_DIR}/H2B5Kac_R1.bigWig \
 ${WWW_DIR}/H2B5Kac_R2.bigWig \
 --labels H2B5Kac_CutNTag_1 H2B5Kac_CutNTag_2 \
 -out H2B5Kac_scores_per_bin.npz \
 --outRawCounts H2B5Kac_scores_per_bin.tab

plotCorrelation \
-in H2B5Kac_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o H2B5Kac_scatter_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix H2B5Kac_scatter_PearsonCorr_bigwigScores.tab

bigwigCompare --operation mean --bigwig1 ${WWW_DIR}/H2B5Kac_R1.bigWig --bigwig2 ${WWW_DIR}/H2B5Kac_R2.bigWig --numberOfProcessors 40 --outFileName ${WWW_DIR}/H2B5Kac.bigWig --outFileFormat bigwig
```

### 2c. Compare all lysine modifications together

#### With all replicates

```console
export WWW_DIR="$PROCESSED/A485/integration/genome_browser/trackHub/mm10/"

multiBigwigSummary bins -b ${WWW_DIR}/As4.1_H3K27Ac_DMSO_REP1.bigWig \
 ${WWW_DIR}/As4.1_H3K27Ac_DMSO_REP2.bigWig \
 ${WWW_DIR}/H3K27ac_R1.bigWig \
 ${WWW_DIR}/H3K27ac_R2.bigWig \
 ${WWW_DIR}/H3K4me1_R1.bigWig \
 ${WWW_DIR}/H3K4me1_R2.bigWig \
 ${WWW_DIR}/H3K4me3_R1.bigWig \
 ${WWW_DIR}/H3K4me3_R2.bigWig \
 ${WWW_DIR}/H4K16ac_R1.bigWig \
 ${WWW_DIR}/H4K16ac_R2.bigWig \
 ${WWW_DIR}/H4K5ac_R1.bigWig \
 ${WWW_DIR}/H4K5ac_R2.bigWig \
 ${WWW_DIR}/H2B5Kac_R1.bigWig \
 ${WWW_DIR}/H2B5Kac_R2.bigWig \
 ${WWW_DIR}/As4.1_ChIP_pooled.bigWig \
 --labels H3K27Ac_DMSO_1 H3K27Ac_DMSO_2 \
 H3K27Ac_CutNTag_1 H3K27Ac_CutNTag_2 \
 H3K4me1_CutNTag_1 H3K4me1_CutNTag_2 \
 H3K4me3_CutNTag_1 H3K4me3_CutNTag_2 \
 H4K16ac_CutNTag_1 H4K16ac_CutNTag_2 \
 H4K5ac_CutNTag_1 H4K5ac_CutNTag_2 \
 H2B5Kac_CutNTag_1 H2B5Kac_CutNTag_2 \
 input \
 -out LysineModifications_scores_per_bin.npz \
 --outRawCounts LysineModifications_scores_per_bin.tab

# Spearman
plotCorrelation \
-in LysineModifications_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o LysineModifications_scatter_SpearmanCorr_bigwigScores.svg   \
--outFileCorMatrix LysineModifications_scatter_SpearmanCorr_bigwigScores.tab

plotCorrelation \
-in LysineModifications_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Average Scores Per Bin" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o LysineModifications_heatmap_SpearmanCorr_bigwigScores.svg \
--outFileCorMatrix LysineModifications_heatmap_SpearmanCorr_bigwigScores.tab

# Pearson 
plotCorrelation \
-in LysineModifications_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--removeOutliers \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o LysineModifications_scatter_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix LysineModifications_scatter_PearsonCorr_bigwigScores.tab

plotCorrelation \
-in LysineModifications_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--removeOutliers \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o LysineModifications_heatmap_PearsonCorr_bigwigScores.svg \
--outFileCorMatrix LysineModifications_heatmap_PearsonCorr_bigwigScores.tab
```

#### Compare after merging replicates

```console
export WWW_DIR="$PROCESSED/A485/integration/genome_browser/trackHub/mm10/"

multiBigwigSummary bins -b ${WWW_DIR}/H3K27Ac_As4.1.bigWig \
 ${WWW_DIR}/H3K4me1.bigWig \
 ${WWW_DIR}/H3K4me3.bigWig \
 ${WWW_DIR}/H4K16ac.bigWig \
 ${WWW_DIR}/H4K5ac.bigWig \
 ${WWW_DIR}/H2B5Kac.bigWig \
 ${WWW_DIR}/As4.1_ChIP_pooled.bigWig \
 --labels H3K27Ac H3K4me1 H3K4me3 \
 H4K16ac H4K5ac H2B5Kac input \
 -out LysineModifications_merged_scores_per_bin.npz \
 --outRawCounts LysineModifications_merged_scores_per_bin.tab

# Spearman
plotCorrelation \
-in LysineModifications_merged_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o LysineModifications_merged_scatter_SpearmanCorr_bigwigScores.svg   \
--outFileCorMatrix LysineModifications_merged_scatter_SpearmanCorr_bigwigScores.tab

plotCorrelation \
-in LysineModifications_merged_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Average Scores Per Bin" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o LysineModifications_merged_heatmap_SpearmanCorr_bigwigScores.svg \
--outFileCorMatrix LysineModifications_merged_heatmap_SpearmanCorr_bigwigScores.tab

# Pearson 
plotCorrelation \
-in LysineModifications_merged_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--removeOutliers \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o LysineModifications_merged_scatter_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix LysineModifications_merged_scatter_PearsonCorr_bigwigScores.tab

plotCorrelation \
-in LysineModifications_merged_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--removeOutliers \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o LysineModifications_merged_heatmap_PearsonCorr_bigwigScores.svg \
--outFileCorMatrix LysineModifications_merged_heatmap_PearsonCorr_bigwigScores.tab

# No input
multiBigwigSummary bins -b ${WWW_DIR}/H3K27Ac_As4.1.bigWig \
 ${WWW_DIR}/H3K4me1.bigWig \
 ${WWW_DIR}/H3K4me3.bigWig \
 ${WWW_DIR}/H4K16ac.bigWig \
 ${WWW_DIR}/H4K5ac.bigWig \
 ${WWW_DIR}/H2B5Kac.bigWig \
 --labels H3K27Ac H3K4me1 H3K4me3 \
 H4K16ac H4K5ac H2B5Kac \
 -out LysineModifications_merged-no-input_scores_per_bin.npz \
 --outRawCounts LysineModifications_merged-no-input_scores_per_bin.tab
 
# Spearman
plotCorrelation \
-in LysineModifications_merged-no-input_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o LysineModifications_merged-no-input_scatter_SpearmanCorr_bigwigScores.svg   \
--outFileCorMatrix LysineModifications_merged-no-input_scatter_SpearmanCorr_bigwigScores.tab

plotCorrelation \
-in LysineModifications_merged-no-input_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Average Scores Per Bin" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o LysineModifications_merged-no-input_heatmap_SpearmanCorr_bigwigScores.svg \
--outFileCorMatrix LysineModifications_merged-no-input_heatmap_SpearmanCorr_bigwigScores.tab

# Pearson 
plotCorrelation \
-in LysineModifications_merged-no-input_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--removeOutliers \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot scatterplot \
-o LysineModifications_merged-no-input_scatter_PearsonCorr_bigwigScores.svg   \
--outFileCorMatrix LysineModifications_merged-no-input_scatter_PearsonCorr_bigwigScores.tab

plotCorrelation \
-in LysineModifications_merged-no-input_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--removeOutliers \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o LysineModifications_merged-no-input_heatmap_PearsonCorr_bigwigScores.svg \
--outFileCorMatrix LysineModifications_merged-no-input_heatmap_PearsonCorr_bigwigScores.tab
```

### 2d. Visualize lysine modifications after scaling

 - mm10 effective genome size: 2,730,871,774

```console
cd $PROCESSED/lysine_modifications/02_alignment/bowtie2/target/markdup/

samtools merge -@ 40 -s 99 --write-index H4K16ac.bam H4K16ac_R1.target.markdup.sorted.bam H4K16ac_R2.target.markdup.sorted.bam

samtools merge -@ 40 -s 99 --write-index H4K5ac.bam H4K5ac_R1.target.markdup.sorted.bam H4K5ac_R2.target.markdup.sorted.bam

samtools merge -@ 40 -s 99 --write-index H3K4me1.bam H3K4me1_R1.target.markdup.sorted.bam H3K4me1_R2.target.markdup.sorted.bam

samtools merge -@ 40 -s 99 --write-index H3K4me3.bam H3K4me3_R1.target.markdup.sorted.bam H3K4me3_R2.target.markdup.sorted.bam

samtools merge -@ 40 -s 99 --write-index H2B5Kac.bam H2B5Kac_R1.target.markdup.sorted.bam H2B5Kac_R2.target.markdup.sorted.bam

samtools merge -@ 40 -s 99 --write-index H3K27ac.bam H3K27ac_R1.target.markdup.sorted.bam H3K27ac_R2.target.markdup.sorted.bam

samtools view -F 260 H4K16ac.bam

bamCoverage --bam H4K16ac.bam -o H4K16ac_RPGC.bw \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2730871774 \
    --ignoreForNormalization chrX \
    --extendReads \
    --numberOfProcessors 40
    
bamCoverage --bam H4K5ac.bam -o H4K5ac_RPGC.bw \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2730871774 \
    --ignoreForNormalization chrX \
    --extendReads \
    --numberOfProcessors 40

bamCoverage --bam H3K4me1.bam -o H3K4me1_RPGC.bw \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2730871774 \
    --ignoreForNormalization chrX \
    --extendReads \
    --numberOfProcessors 40

bamCoverage --bam H3K4me3.bam -o H3K4me3_RPGC.bw \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2730871774 \
    --ignoreForNormalization chrX \
    --extendReads \
    --numberOfProcessors 40

bamCoverage --bam H2B5Kac.bam -o H2B5Kac_RPGC.bw \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2730871774 \
    --ignoreForNormalization chrX \
    --extendReads \
    --numberOfProcessors 40
    
bamCoverage --bam H3K27ac.bam -o H3K27ac_RPGC.bw \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2730871774 \
    --ignoreForNormalization chrX \
    --extendReads \
    --numberOfProcessors 40

cp *_RPGC.bw $PROCESSED/A485/integration/genome_browser/trackHub/mm10/
```

### 2e. Move files to public genome_browser/ locations

```console
export PROJECT_WWW="$WWW/genome_browser/trackHub/mm10/"
export WWW_DIR="$PROCESSED/A485/integration/genome_browser/trackHub/mm10/"

cp ${PROJECT_WWW}/H3K27ac_consensus_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/H4K16ac_consensus_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/H4K5ac_consensus_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/H2B5Kac_consensus_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/H3K4me3_consensus_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/H3K4me1_consensus_peaks.bb ${WWW_DIR}
cp ${PROJECT_WWW}/H2B5Kac_R1.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H2B5Kac_R2.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H4K5ac_R1.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H4K5ac_R2.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H3K4me3_R1.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H3K4me3_R2.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H3K4me1_R1.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H3K4me1_R2.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H3K27ac_R1.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H3K27ac_R2.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H4K16ac_R1.bigWig ${WWW_DIR}
cp ${PROJECT_WWW}/H4K16ac_R2.bigWig ${WWW_DIR}
```