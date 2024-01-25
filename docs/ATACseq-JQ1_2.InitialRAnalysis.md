
## 2. Load R and begin primary analysis

```R
# devtools::install_github("andrelmartins/bigwig")
library(DESeq2)
library(lattice)
library(bigWig)
library(genefilter)
library(ggplot2)
library(ggrepel)
library(svglite)
library(DEGreport)
library(tibble)
```

### 2a. Set up helper functions

```R
#functions
getRawCountsInterval <- function(df, bigWig_path, file_prefix = 'M') {
	df = df[,1:3]
	vec_names <- c()
	inten_df  <- data.frame(matrix(ncol = 0, nrow = nrow(df)))
	for (mod_bigWig in Sys.glob(file.path(bigWig_path,
            paste(file_prefix, "*.bigWig", sep ='')))) {
		factor_name <- strsplit(strsplit(mod_bigWig, "/")[[1]][length(strsplit(mod_bigWig, "/")[[1]])], '\\.')[[1]][1]
		print(factor_name)
		vec_names <- c(vec_names, factor_name)
		loaded_bw <- load.bigWig(mod_bigWig)
		mod_inten <- bed.region.bpQuery.bigWig(loaded_bw, df)
		inten_df  <- cbind(inten_df, mod_inten)
	}
	colnames(inten_df) <- vec_names
	rowname_vals <- paste(df[,1], ':', df[,2], '-', df[,3], sep='')
	row.names(inten_df) <- rowname_vals
	return(inten_df)
}
```

### 2b. Perform initial processing of raw ATAC data

```R
#processes
directory = "$PROCESSED/jq1/bam_files/"
setwd(directory)

#Change the file to the one you want to read
consensus_peaks_file <- read.table("$PROCESSED/jq1/bam_files/As4.1_smallpval_rmBlacklist.bed")

bigWig_dir <- "$PROCESSED/jq1/bam_files/submission"
trt_JQ1_counts   <- getRawCountsInterval(consensus_peaks_file,
                                         bigWig_dir, file_prefix = 'JQ1')
ctrl_DMSO_counts <- getRawCountsInterval(consensus_peaks_file,
                                         bigWig_dir, file_prefix = 'DMSO')

counts_df <- cbind(trt_JQ1_counts, ctrl_DMSO_counts)
save(counts_df, file='$PROCESSED/jq1/bam_files/counts_df.Rdata')

sample_conditions = c("JQ1", "JQ1", "DMSO", "DMSO")
deseq_counts_table = DESeqDataSetFromMatrix(counts_df,
    DataFrame(sample_conditions), ~ sample_conditions)
colData(deseq_counts_table)$condition <- c("JQ1", "JQ1", "DMSO", "DMSO")
dds = DESeq(deseq_counts_table)
counts_df_deseq_raw  <- dds

#counts table
normalized_counts_atac = counts(counts_df_deseq_raw, normalized=TRUE)
save(normalized_counts_atac, file='normalized_counts_atac.Rdata')

#PCA
rld <- rlog(counts_df_deseq_raw, blind=TRUE)
save(rld, file = '$PROCESSED/jq1/bam_files/rld.Rdata')
#clustering
ddsLRT <- DESeq(counts_df_deseq_raw, test="LRT", reduced = ~ 1)

resLRT <- results(ddsLRT)
save(resLRT, file = 'resLRT.Rdata')

resLRT[!is.na(resLRT$padj),]
# DataFrame with 106078 rows and 6 columns
range(resLRT[!is.na(resLRT$padj),]$padj)
summary(resLRT[!is.na(resLRT$padj),]$padj)

padj.cutoff <- 0.00000001 #1e-8 362 rows
padj.cutoff <- 0.000001 #1e-6 726 rows
padj.cutoff <- 0.00001 #1e-5 1093 rows
padj.cutoff <- 0.001 #1e-3 3195 rows

siglrtRE <- resLRT[resLRT$padj < padj.cutoff & !is.na(resLRT$padj),]

rld_mat        <- assay(rld)
cluster_rlog   <- rld_mat[rownames(siglrtRE),]

meta           <- data.frame(sample_conditions=c("JQ1_01", "JQ1_02", "DMSO_01", "DMSO_02"))
rownames(meta) <- colnames(cluster_rlog)
meta$sample_conditions <- factor(meta$sample_conditions,
    levels=c("JQ1_01", "JQ1_02", "DMSO_01", "DMSO_02"))
meta$group <- factor(c("JQ1", "JQ1", "DMSO", "DMSO"))
save(cluster_rlog, meta, sample_conditions, file = 'cluster_rlog_pval_1e8.Rdata')
```

```R
#generate nondynamic peaks set for FIMO
not_different <- rownames(resLRT[resLRT$padj > 0.5 & !is.na(resLRT$padj) & !is.na(resLRT$log2FoldChange) & abs(resLRT$log2FoldChange) < 0.25,])
chr   <- sapply(strsplit(not_different, ':'), '[', 1)
x     <- sapply(strsplit(not_different, ':'), '[', 2)
start <- sapply(strsplit(x, '-'), '[', 1)
end   <- sapply(strsplit(x, '-'), '[', 2)

curated_not_different <- data.frame(chr,start,end)
write.table(curated_not_different, file='nondynamic_peaks.bed',
            sep='\t', col.names=F, row.names=F, quote=F)

#generate all peaks set
chr   <- sapply(strsplit(rownames(resLRT), ':'), '[', 1)
z     <- sapply(strsplit(rownames(resLRT), ':'), '[', 2)
start <- sapply(strsplit(z, '-'), '[', 1)
end   <- sapply(strsplit(z, '-'), '[', 2)

bed <- data.frame(chr, start ,end, resLRT$baseMean)
write.table(bed[,1:3], file = paste0('all_peaks_', Sys.Date(), '.bed'),
            sep='\t', col.names=F, row.names=F, quote=F)
```

### 2c. Plot nondynamic vs dynamic peaks

```R
plot_df <- data.frame(type=c("Dynamic", "Nondynamic"),
                      num=c(nrow(cluster_rlog),
                            (nrow(resLRT)-nrow(cluster_rlog))),
                      x=1)
svglite(paste0("num_dynamic-vs-nondynamic_peaks_", Sys.Date(), ".svg"), 
        bg = "transparent", fix_text_size = FALSE, pointsize = 6)
pdf(paste0("num_dynamic-vs-nondynamic_peaks_", Sys.Date(), ".pdf"))      
ggplot(plot_df, aes(x=x, y=num, fill=type)) +
    geom_bar(stat='Identity', color='black') +
    labs(title = NULL,
    y = 'Number of Peaks',
    x = NULL,
    fill = 'Type') +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=16,face='bold',color='black'),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=18,face='bold'),
          legend.title = element_text(size=16,face='bold'),
          legend.text = element_text(size=14,face='bold'),
          plot.title = element_text(size=17,face='bold',hjust=0.5)) +
    scale_fill_manual(values = c('#f16469','#f1bc7b'))
dev.off()
```

### 2d. Perform dynamic peak clustering

```R
svglite(paste0("degPatterns_1e3_", Sys.Date(), ".svg"),
        bg = "transparent", fix_text_size = FALSE, pointsize = 6)
pdf(paste0("degPatterns_1e3_", Sys.Date(), ".svg"))
clusters_all_1e3 <- DEGreport::degPatterns(cluster_rlog, metadata = meta,
                                minc = 100, time = "sample_conditions",
                                col=NULL, eachStep = FALSE)
dev.off()
save(clusters_all_1e3, file = paste0('clusters_all_minc100_1e3_', Sys.Date(), '.Rdata')
```
Cluster1: INCREASED accessibility with JQ1 treatment (1417 regions)
Cluster2: DECREASED accessibility with JQ1 treatment (1778 regions)

Test higher p cutoff
```R
padj.cutoff <- 0.00000001 #1e-8 362 rows
padj.cutoff <- 0.000001 #1e-6 726 rows
padj.cutoff_1e5 <- 0.00001 #1e-5 1093 rows
siglrtRE_e5 <- resLRT[resLRT$padj < padj.cutoff_1e5 & !is.na(resLRT$padj),]
rld_mat_e5        <- assay(rld)
cluster_rlog_e5   <- rld_mat_e5[rownames(siglrtRE_e5),]
svglite(paste0("degPatterns_1e5_", Sys.Date(), ".svg"),
        bg = "transparent", fix_text_size = FALSE, pointsize = 6)
clusters_all_1e5 <- DEGreport::degPatterns(cluster_rlog_e5, metadata = meta,
                                minc = 100, time = "sample_conditions",
                                col=NULL, eachStep = FALSE)
dev.off()
save(clusters_all_1e5, file = paste0('clusters_all_minc100_1e5_', Sys.Date(), '.Rdata'))
```
1e-05
Cluster1: INCREASED accessibility with JQ1 treatment (468 regions)
Cluster2: DECREASED accessibility with JQ1 treatment (625 regions)

Plotting
```R
#generate 'plot.df' object
plot_df = clusters_all_1e3$normalized
plot_df$cutoff0.444[is.na(plot_df$cutoff0.444)] <- 0

plot_df$sample_conditions <- as.character(plot_df$sample_conditions)
plot_df$sample_conditions[plot_df$sample_conditions == 'DMSO_01'] <- 0
plot_df$sample_conditions[plot_df$sample_conditions == 'DMSO_02'] <- 1
plot_df$sample_conditions[plot_df$sample_conditions == 'JQ1_01']  <- 2
plot_df$sample_conditions[plot_df$sample_conditions == 'JQ1_02']  <- 3
plot_df$sample_conditions <- as.numeric(plot_df$sample_conditions)
plot_df <- plot_df[order(plot_df$genes),]
plot_df <- plot_df[order(plot_df$sample_conditions),]

plot_df$cluster <- paste('cluster', as.character(plot_df$cluster), sep = '')

plot_df$chr   <- sapply(strsplit(plot_df$genes, '[.]'), '[', 1)
plot_df$start <- sapply(strsplit(plot_df$genes, '[.]'), '[', 2)
plot_df$end   <- sapply(strsplit(plot_df$genes, '[.]'), '[', 3)

write.table(cbind(plot_df$chr, plot_df$start, plot_df$end),
            file = paste0('dynamic_peaks_', Sys.Date(), '.bed'),
            quote = FALSE, sep = '\t', col.names=FALSE, row.names=FALSE)
```

### 2e. Save a BED file for each unique cluster.

```R
for (i in unique(plot_df$cluster)) {
    print(i)
    write.table(plot_df[plot_df$cluster == i,
                        c('chr','start','end', 'value', 'cluster')][!duplicated(plot_df[plot_df$cluster == i,]$genes),],
                file = paste0('cluster_bed_',
                              gsub(" ", "", i, fixed = TRUE), '_',
                              Sys.Date(), '.bed'),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
}

plot_df$supercluster <- 'down'
plot_df[plot_df$cluster %in% c("cluster2"),'supercluster'] <- 'up'

#number of peaks in each supercluster
table(plot_df$supercluster) / 2
plot_df_atac = plot_df[,c(1,3,4,6,11:14)]
plot_df_atac$genes = paste0(plot_df_atac$chr, ':', plot_df_atac$start, '-', plot_df_atac$end)
colnames(plot_df_atac)[1] <- 'peak'
colnames(plot_df_atac)[3] <- 'treatment'
save(plot_df_atac, file='plot_df_atac.Rdata')
```