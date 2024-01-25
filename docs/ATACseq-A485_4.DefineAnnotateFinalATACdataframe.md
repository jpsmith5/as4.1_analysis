### 4a. Create ATAC data frame

```R
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(DESeq2)
library(lattice)
library(bigWig)
library(ggplot2)
library(ggrepel)
library(svglite)
library(DEGreport)
library(tibble)
library(fastSave)

directory = "$PROCESSED/A485/bam_files/"
tomtom_dir = "$PROCESSED/A485/bam_files/FIMO/tomtom/"

# counts table
head(normalized_counts_atac)
categorize_deseq_df <- function(df, fdr = 0.05, log2fold = 0.0, treat = 'A485') {
    df_effects_lattice = df
    df_effects_lattice$response = 'All Other Peaks'
    if(nrow(df_effects_lattice[df_effects_lattice$padj < fdr & !is.na(df_effects_lattice$padj) & df_effects_lattice$log2FoldChange > log2fold,]) > 0) {
        df_effects_lattice[df_effects_lattice$padj < fdr & !is.na(df_effects_lattice$padj) & df_effects_lattice$log2FoldChange > log2fold,]$response = 'Increased'
    }
    if(nrow(df_effects_lattice[df_effects_lattice$padj < fdr & !is.na(df_effects_lattice$padj) & df_effects_lattice$log2FoldChange < log2fold,]) > 0) {
        df_effects_lattice[df_effects_lattice$padj < fdr & !is.na(df_effects_lattice$padj) & df_effects_lattice$log2FoldChange < log2fold,]$response = 'Decreased'
    }
    if(nrow(df_effects_lattice[df_effects_lattice$padj > 0.5 & !is.na(df_effects_lattice$padj) & abs(df_effects_lattice$log2FoldChange) < 0.25,]) > 0) {
        df_effects_lattice[df_effects_lattice$padj > 0.5 & !is.na(df_effects_lattice$padj) & abs(df_effects_lattice$log2FoldChange) < 0.25,]$response = 'Unchanged'
    }
    return(df_effects_lattice)
}

peaks = unique(read.table(file.path(directory, "all_peaks.bed"),sep='\t'))
peaks$location = paste0(peaks$V1,':',peaks$V2,'-',peaks$V3)
df = data.frame(row.names = peaks$location, chr=peaks[,1], start=peaks[,2], end=peaks[,3])

control = c()
treated = c()
genes = c()

print('Getting ATAC means')
for (i in 1:nrow(normalized_counts_atac)) {
    control = append(control, mean(normalized_counts_atac[i,3:4])) # DMSO_01, DMSO_02
    treated = append(treated, mean(normalized_counts_atac[i,1:2])) # A485_01, A485_02
    genes = append(genes,rownames(normalized_counts_atac)[i])
}

atac_means = data.frame(row.names = genes, ctrl = control, trt = treated)
head(atac_means)

df_merged = merge(df, atac_means, by='row.names', all=TRUE)
rownames(df_merged) = df_merged[,1]
df_merged = df_merged[,-1]
head(df_merged)
save(df_merged, file='df_merged_10.23.2023.RData')

load(file=file.path(tomtom_dir, 'df_merged.RData'))

cluster1_bed <- fread(file.path(directory, "cluster_bed_cluster1.bed"))
# cluster 1 is the group that goes UP with treatment (A485)

cluster2_bed <- fread(file.path(directory, "cluster_bed_cluster2.bed"))
# cluster 2 is the group that goes DOWN with treatment (A485)
```

### 4b. Modify BED files for use in GREAT.

```R
cluster1_bed_reformat <- data.table(
    chr=cluster1_bed$V1,
    start=cluster1_bed$V2,
    end=cluster1_bed$V3,
    name=paste0("cluster1_", seq.int(1, nrow(cluster1_bed), 1)),
    score=round(scales::rescale(cluster1_bed$V4, c(0, 1000))),
    strand= ".")
rtracklayer::export.bed(cluster1_bed_reformat, "cluster_bed_cluster1.bed6")

cluster2_bed_reformat <- data.table(
    chr=cluster2_bed$V1,
    start=cluster2_bed$V2,
    end=cluster2_bed$V3,
    name=paste0("cluster2_", seq.int(1, nrow(cluster2_bed), 1)),
    score=round(scales::rescale(cluster2_bed$V4, c(0, 1000))),
    strand= ".")
rtracklayer::export.bed(cluster2_bed_reformat, "cluster_bed_cluster2.bed6")

plot_dt <- as.data.table(plot_df)
plot_dt <- splitstackshape::getanID(data = plot_dt, id.vars = "cluster")

# this is all dynamic peaks, which includes duplicate regions (same set of regions for each input sample, A485_01, 02, and DMSO_01, 02.

# drop the duplicates
plot_dt <- plot_dt[!duplicated(plot_df$genes),]
dynamic_peaks_reformat <- data.table(
    chr=plot_dt$chr,
    start=plot_dt$start,
    end=plot_dt$end,
    name=paste0(plot_dt$cluster, "_", plot_dt$.id),
    score=round(scales::rescale(plot_dt$value, c(0, 1000))),
    strand= ".")

# Or, just combined cluster1 and cluster2
dynamic_peaks_background <- rbind(cluster1_bed_reformat, cluster2_bed_reformat)
dynamic_peaks_background <- dynamic_peaks_background[order(chr, start, end),]
rtracklayer::export.bed(dynamic_peaks_background, "cluster_bed_background.bed6")
```

```R
#add cluster information
print('Adding Cluster Info')
print(head(df_merged))

supercluster_df = data.frame(cluster=c(1,2),
                             supercluster=c(rep('up',1),
                                            rep('down',1)))
# THis step needs more memory, come back to it.
cluster_df = data.frame()
for (file in Sys.glob(file.path(directory, 'cluster_bed_cluster*.bed'))) {
    cluster = as.numeric(strsplit(strsplit(file,'cluster_bed_cluster')[[1]][2],'.bed')[[1]][1])
    print(cluster)
    sc = supercluster_df[supercluster_df$cluster == cluster,]$supercluster
    x = read.table(file)
    x$location = paste0(x$V1,':',x$V2,'-',x$V3)
    cluster_df = rbind(cluster_df, data.frame(peak=x$location, 
                       cluster=cluster, supercluster=sc))
}

# my method is faster
dt1 <- as.data.table(df_merged, keep.rownames = TRUE, key = 'rn')
dt1 <- dt1[order(chr, start, end),]
dt2 <- as.data.table(cluster_df)
dt1_subset <- dt1[rn %in% dt2$peak,]
dt_merged_cluster <- merge(dt1_subset, dt2, by.x="rn", by.y="peak", all.x=TRUE)

# arun's method
df_merged_cluster = merge(df_merged, cluster_df,
                          by.x='row.names', by.y=1, all.x=TRUE)

rownames(df_merged_cluster) = df_merged_cluster[,1]
df_merged_cluster = df_merged_cluster[,-1]

save(df_merged_cluster, file='df_merged_cluster.RData')
```

```R
#add tf family scores
print('Adding TF Family Scores')
print(head(df_merged_cluster))

setwd("$PROCESSED/A485/bam_files/FIMO/tomtom/")

scores_df = data.frame()
for (file in Sys.glob('fimo_composites/main_figure_beds/*_fimo_all.bed')) {
    factor = strsplit(strsplit(file,'/')[[1]][3],'_fimo_all.bed')[[1]][1]
    print(factor)
    x = read.table(file,sep='\t')
    x = x[x$V5 != -1,]
    if (nrow(x) !=0) {
        if(factor %in% c('Sp','Klf')) {
            y = aggregate(as.numeric(V7)~V1+V2+V3, data=x, FUN=sum)
        } else {
            y = aggregate(as.numeric(V8)~V1+V2+V3, data=x, FUN=sum)
        }
        colnames(y) = c('chr', 'start', 'end', factor)
        rownames(y) = paste0(y[,1], ':', y[,2], '-', y[,3])
        
        scores_df = rbind(scores_df, data.frame(peak = rownames(y),
                          score=y[,4], factor = factor))
    } else {
        paste0("WARNING: ", factor, " does not contain positive scores.")
    }
}

table(scores_df$factor)
```

 Fra2  Klf4   Maz  RUNX   Sp2 TEAD4 Trp63
59171 37604 21229 28958 38557 27021 20963

```R
fimo_scores_all_atac = data.frame(peak = unique(scores_df$peak))

for(factor in unique(scores_df$factor)) {
    temp = scores_df[scores_df$factor == factor,]
    temp = temp[,c(1,2)]
    fimo_scores_all_atac = merge(fimo_scores_all_atac,temp,by='peak',all.x=TRUE)
    colnames(fimo_scores_all_atac)[ncol(fimo_scores_all_atac)] = factor
}

rownames(fimo_scores_all_atac) = fimo_scores_all_atac$peak
fimo_scores_all_atac = fimo_scores_all_atac[,-1]
save(fimo_scores_all_atac, file='fimo_scores_all_atac.RData')

load('fimo_scores_all_atac.RData')

fimo_scores_all_atac$peak <- rownames(fimo_scores_all_atac)
df_scores <- merge(dt_merged_cluster, fimo_scores_all_atac,
                   by.x="rn", by.y="peak", all.x=TRUE)
                  
df_scores = df_scores[,-1]
rownames(df_scores) = dt_merged_cluster$rn
save(df_scores, file='df_scores.RData')
```

```R
#add pairwise comparisons
print('Adding Pairwise Comparisons')
print(head(df_scores))

condition_key = data.frame(condition=c(rep('A485', 2),
                                       rep('DMSO', 2)),
                                       cols = c(1:4))
condition = c('A485', 'DMSO')

comparisons_df = data.frame()
i <- 1
j <- 2

print(paste0(condition[j],'.v.',condition[i]))

a = condition_key[condition_key$condition == condition[i],]$cols
b = condition_key[condition_key$condition == condition[j],]$cols
merged_counts_small = counts_df[,c(a,b)]

# number of replicates per condition
unt = 2
trt = 2

sample_conditions = factor(c(rep("treated",trt),
                             rep("untreated",unt)),
                           levels=c("untreated", "treated"))
mm_deseq_counts_table = DESeqDataSetFromMatrix(merged_counts_small,
    DataFrame(sample_conditions), ~ sample_conditions)

mm_atac = mm_deseq_counts_table
atac_size_factors = estimateSizeFactorsForMatrix(merged_counts_small)

sizeFactors(mm_atac) = atac_size_factors
mm_atac = estimateDispersions(mm_atac)
mm_atac = nbinomWaldTest(mm_atac)
res_mm_atac = results(mm_atac)

lattice = categorize_deseq_df(res_mm_atac,
                              fdr = 0.001,
                              log2fold = 0.0,
                              treat = '')
table(lattice$response)

<!-- 
All Other Peaks       Decreased       Increased       Unchanged
          89372            4771            4307           29796
-->
lattice = as.data.frame(lattice[,c(2,6,7)])
colnames(lattice) = paste0(colnames(lattice), '.',
                           condition[j], '.v.', condition[i])

comparisons_df = merge(comparisons_df, lattice, by='row.names', all=TRUE)
rownames(comparisons_df) = comparisons_df$Row.names
comparisons_df = comparisons_df[,-1]

save(comparisons_df, file='comparisons_df.RData')

comparisons_df$peak <- rownames(comparisons_df)
df_scores$rn <- rownames(df_scores)
df_scores_comparisons <- merge(df_scores, comparisons_df,
                               by.x="rn", by.y="peak", all.x=TRUE)
df_scores_comparisons <- as.data.frame(df_scores_comparisons)
rownames(df_scores_comparisons) = df_scores_comparisons$rn
head(df_scores_comparisons)

df_scores_comparisons = df_scores_comparisons[,-1]
head(df_scores_comparisons)
save(df_scores_comparisons, file='df_scores_comparisons.RData')
```

### 4c. Annotate genomic features

```R
print('Add baseMean and overall time course info')

padj_cutoff <- 0.00000001
padj_cutoff <- 0.00001 #e-05

load('resLRT.RData')
res_lrt = as.data.frame(resLRT[,c(1,2,6)])
res_lrt$response = 'Nondynamic'
res_lrt[res_lrt$padj < padj_cutoff & !is.na(res_lrt$padj),]$response = 'Dynamic'
colnames(res_lrt) = paste0(colnames(res_lrt),'_DMSOvA485')
res_lrt$peak <- rownames(res_lrt)

#df_lrt = merge(df_scores_comparisons, res_lrt, by = 'row.names', all = TRUE)
df_scores_comparisons$rn <- rownames(df_scores_comparisons)
df_lrt <- merge(df_scores_comparisons, res_lrt, by.x="rn", by.y="peak", all.x=TRUE)
head(df_lrt)

rownames(df_lrt) = paste0(df_lrt$chr, ':', df_lrt$start, '-', df_lrt$end)
df_lrt = df_lrt[,-1]

#add distribution
print('Add distribution')

df_lrt$distribution = 'Intergenic'

library(org.Mm.eg.db)
library(ChIPpeakAnno)
data(TSS.mouse.NCBIM37)
data(TSS.mouse.GRCm38)

all_peaks <- fread(file.path(directory, 'all_peaks.bed'))
colnames(all_peaks) <- c("chr", "start", "end")
all_peaks_GR <- makeGRangesFromDataFrame(all_peaks, keep.extra.columns=T)
annotated_peaks <- ChIPpeakAnno::annotatePeakInBatch(all_peaks_GR, 
    AnnotationData = TSS.mouse.GRCm38)
annotated_peaks <- ChIPpeakAnno::addGeneIDs(annotated_peaks,
    "org.Mm.eg.db",
    c("symbol","genename"))
annotated_peaks_dt <- as.data.table(annotated_peaks)
fwrite(annotated_peaks_dt, "As4.1_A485-vs-DMSO_all_peaks_annotated.csv", sep=",")

library(reactome.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Peak distribution over genomic features
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
all_peaks_featuresDist <- assignChromosomeRegion(annotated_peaks,
    nucleotideLevel=FALSE,
    precedence=c("Promoters", "immediateDownstream", "fiveUTRs",
                 "threeUTRs","Exons", "Introns"), TxDb=txdb)

svglite("all_peaks_featuresDist_A485_vs_DMSO.svg")
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(all_peaks_featuresDist$percentage, las=1, horiz=T)
dev.off()

annotated_peaks2 <- genomicElementDistribution(annotated_peaks, TxDb=txdb, nucleotideLevel = FALSE, plot=FALSE)
annotated_peaks_dt2 <- as.data.table(annotated_peaks2$peaks)
fwrite(annotated_peaks_dt2, "As4.1_A485-vs-DMSO_all_peaks_annotated_genomic_region.csv", sep=",")

svglite("all_peaks_pie_A485_vs_DMSO.svg")
genomicElementDistribution(annotated_peaks, TxDb=txdb, nucleotideLevel = FALSE, plot=TRUE)
dev.off()

library(UpSetR)
upset_dt <- genomicElementUpSetR(annotated_peaks, TxDb=txdb)

svglite("all_peaks_upset_A485_vs_DMSO.svg", width=20, 
        bg = "transparent", fix_text_size = FALSE, pointsize = 6)
upset(upset_dt$plotData, nsets=13, nintersects=NA)
dev.off()

annotated_peaks_dt2[, region := paste0(seqnames, ":", start, "-", end)]
dt_lrt <- as.data.table(df_lrt, keep.rownames=T, key="rn")
annotated_dt <- merge(dt_lrt, annotated_peaks_dt2, by.x="rn", by.y="region")
table(annotated_dt$geneLevel, annotated_dt$ExonIntron)

<!--
                   exon intergenic intron
  distalIntergenic    0       1836      0
  geneBody          231         31   1597
  geneDownstream      2         42      7
  promoter          106         95     23
-->

final_atac_all_dt <- annotated_dt
fwrite(final_atac_all_dt, paste0("final_atac_annotated_data_", Sys.Date(), ".csv"), sep=",")
save(final_atac_all_dt, file=paste0('final_atac_all_dt_', Sys.Date(), '.RData'))
```

Annotate the peaks that make up the UP or DOWN sets (cluster 1 and 2)

```R
annotated_UP <- annotated_dt[response_DMSOvA485=="Dynamic" & supercluster == "up",]

colnames(annotated_UP) <- c("rn", "chr", "start", "end", "ctrl", "trt",
                            "cluster", "supercluster",
                            "Fra2", "Klf4", "Maz", "RUNX", "Sp2", "TEAD4",
                            "Trp63", 
                            "log2FoldChange.DMSO.v.A485", "padj.DMSO.v.A485",
                            "response.DMSO.v.A485", "baseMean_DMSOvA485",
                            "log2FoldChange_DMSOvA485", "padj_DMSOvA485",
                            "response_DMSOvA485", "distribution", "chr2",
                            "start2", "end2", "width", "strand", "peak",
                            "feature", "start_position", "end_position",
                            "feature_strand", "insideFeature",
                            "distancetoFeature", "shortestDistance",
                            "fromOverlappingOrNearest", "symbol",
                            "genename", "geneLevel", "ExonIntron", "Exons")
annotated_UP$`chr2` <- NULL
annotated_UP$`start2` <- NULL
annotated_UP$`end2` <- NULL
annotated_UP_GR <- makeGRangesFromDataFrame(annotated_UP, keep.extra.columns=T)

table(annotated_UP$ExonIntron)
table(annotated_UP$ExonIntron, annotated_UP$geneLevel)
table(annotated_UP$geneLevel)
```

      exon intergenic     intron
       232        843        765

             distalIntergenic geneBody geneDownstream promoter
  exon                      0      140              0       92
  intergenic              738       15             27       63
  intron                    0      742              5       18

distalIntergenic         geneBody   geneDownstream         promoter
             738              897               32              173


```R
annotated_UP_featuresDist <- assignChromosomeRegion(annotated_UP_GR,
    nucleotideLevel=FALSE,
    precedence=c("Promoters", "immediateDownstream", "fiveUTRs",
                 "threeUTRs","Exons", "Introns"), TxDb=txdb)

svglite("cluster1_UP_featuresDist_A485_vs_DMSO.svg")
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(annotated_UP_featuresDist$percentage, las=1, horiz=T)
dev.off()

upset_UP_dt <- genomicElementUpSetR(annotated_UP_GR, TxDb=txdb)

svglite("cluster1_UP_upset_A485_vs_DMSO.svg", width=20, 
        bg = "transparent", fix_text_size = FALSE, pointsize = 6)
upset(upset_UP_dt$plotData, nsets=13, nintersects=NA)
dev.off()

svglite("cluster1_UP_pie_A485_vs_DMSO.svg")
genomicElementDistribution(annotated_UP_GR, TxDb=txdb, nucleotideLevel = FALSE, plot=TRUE)
dev.off()
```

Now, annotate the cluster that goes down with treatment

```R
annotated_DOWN <- annotated_dt[response_DMSOvA485=="Dynamic" & supercluster == "down",]
colnames(annotated_DOWN) <- c("rn", "chr", "start", "end", "ctrl", "trt",
                            "cluster", "supercluster",
                            "Fra2", "Klf4", "Maz", "RUNX", "Sp2", "TEAD4",
                            "Trp63",
                            "log2FoldChange.DMSO.v.A485", "padj.DMSO.v.A485",
                            "response.DMSO.v.A485", "baseMean_DMSOvA485",
                            "log2FoldChange_DMSOvA485", "padj_DMSOvA485",
                            "response_DMSOvA485", "distribution", "chr2",
                            "start2", "end2", "width", "strand", "peak",
                            "feature", "start_position", "end_position",
                            "feature_strand", "insideFeature",
                            "distancetoFeature", "shortestDistance",
                            "fromOverlappingOrNearest", "symbol",
                            "genename", "geneLevel", "ExonIntron", "Exons")
annotated_DOWN$`chr2` <- NULL
annotated_DOWN$`start2` <- NULL
annotated_DOWN$`end2` <- NULL
annotated_DOWN_GR <- makeGRangesFromDataFrame(annotated_DOWN, keep.extra.columns=T)

table(annotated_DOWN$ExonIntron)
table(annotated_DOWN$ExonIntron, annotated_DOWN$geneLevel)
table(annotated_DOWN$geneLevel)
```
> table(annotated_DOWN$ExonIntron)

      exon intergenic     intron
       107       1161        862
> table(annotated_DOWN$ExonIntron, annotated_DOWN$geneLevel)

             distalIntergenic geneBody geneDownstream promoter
  exon                      0       91              2       14
  intergenic             1098       16             15       32
  intron                    0      855              2        5
> table(annotated_DOWN$geneLevel)

distalIntergenic         geneBody   geneDownstream         promoter
            1098              962               19               51

```R
annotated_DOWN_featuresDist <- assignChromosomeRegion(annotated_DOWN_GR,
    nucleotideLevel=FALSE,
    precedence=c("Promoters", "immediateDownstream", "fiveUTRs",
                 "threeUTRs","Exons", "Introns"), TxDb=txdb)

svglite("cluster2_DOWN_featuresDist_A485_vs_DMSO.svg")
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(annotated_DOWN_featuresDist$percentage, las=1, horiz=T)
dev.off()

upset_UP_dt <- genomicElementUpSetR(annotated_DOWN_GR, TxDb=txdb)

svglite("cluster2_DOWN_upset_A485_vs_DMSO.svg", width=20, 
        bg = "transparent", fix_text_size = FALSE, pointsize = 6)
upset(upset_UP_dt$plotData, nsets=13, nintersects=NA)
dev.off()

svglite("cluster2_DOWN_pie_A485_vs_DMSO.svg")
genomicElementDistribution(annotated_DOWN_GR, TxDb=txdb, nucleotideLevel = FALSE, plot=TRUE)
dev.off()
```

Look at the top represented genes based on target regions

```R
top_genes <- as.data.table(sort(table(final_atac_all_dt$symbol)))
colnames(top_genes) <- c("gene", "N")
fwrite(top_genes, "final_atac_gene_counts.csv", sep=",")
```
