###  11. Functional analysis v2 (`clusterProfiler`)

#### 11a. H3K27Ac

```R
project_name <- "H3K27Ac"

loadGRanges <- function(file_name) {
    tmp <- fread(file_name)
    colnames(tmp) <- c("chr", "start", "end", "name", "score", "strand")
    tmp2 <- makeGRangesFromDataFrame(tmp, keep.extra.columns=TRUE)
    seqlevelsStyle(tmp2) <- "UCSC"
    return(tmp2)
}

# Load data
samplefiles <- list.files(paste0(project_name, "_DB/"),
                          pattern= "_true_peaks.bed", full.names=T)
samplefiles <- sapply(file.path(samplefiles), loadGRanges)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("DMSO", "H89")

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)

fwrite(as.data.table(peakAnnoList$DMSO@annoStat),
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_H89_feature_distribution.csv"))

fwrite(as.data.table(peakAnnoList$H89@annoStat),
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_DMSO_feature_distribution.csv"))

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_H89-v-DMSO_feature_distribution_",
               Sys.Date(), ".svg"))
plotAnnoBar(peakAnnoList)
dev.off()

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_H89-v-DMSO_TSS_distribution_",
               Sys.Date(), ".svg"))
plotDistToTSS(peakAnnoList,
    title="Distribution of transcription factor-binding loci \n relative to TSS")
dev.off()

# As a GRanges object, just the peak annotations for H89 for example
H3K27Ac_H89_annot <- data.frame(peakAnnoList[["H89"]]@anno)

# Get the entrez IDs
entrez_IDs <- H3K27Ac_H89_annot$geneId

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                                         keys = entrez_IDs,
                                         columns = c("GENENAME"),
                                         keytype = "ENTREZID")

# Change IDs to character type to merge
annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)

# Write to file
H3K27Ac_H89_annot %>% 
  left_join(annotations_edb, by=c("geneId"="ENTREZID")) %>% 
  write.table(file=paste0(project_name, "_DB/", project_name,
                          "_H89_peak_annotation.txt"),
              sep="\t", quote=F, row.names=F)
  
# Run clusterProfiler's GO enrichment analysis 
ego <- clusterProfiler::enrichGO(gene = entrez_IDs, 
                                 keyType = "ENTREZID", 
                                 OrgDb = org.Mm.eg.db, 
                                 ont = "BP", 
                                 pAdjustMethod = "BH", 
                                 qvalueCutoff = 0.05, 
                                 readable = TRUE)
# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary,
          paste0(project_name, "_DB/", project_name,
                 "_H89_true_peaks_clusterProfiler_GO.csv"))

# Dotplot visualization
svglite(paste0(project_name, "_DB/", project_name,
               "_H89_true_peaks_clusterProfiler_GO_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(ego, showCategory=50) +
    ggplot2::scale_y_discrete(labels=function(x) str_wrap(x, width=70))
dev.off()

# Compare GO analysis
## Create a list with genes from each sample
genes  <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
compGO <- clusterProfiler::compareCluster(geneCluster = genes, 
                                          fun = "enrichGO",
                                          OrgDb = org.Mm.eg.db, 
                                          pvalueCutoff  = 0.05, 
                                          pAdjustMethod = "BH")
svglite(paste0(project_name, "_DB/", project_name,
                "_H89-v-DMSO_true_peaks_clusterProfiler_GO_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(compGO, showCategory = 20,
        title = "GO Pathway Enrichment Analysis") +
    ggplot2::scale_y_discrete(labels=function(x) str_wrap(x, width=70))
dev.off()

# Evaluate and compare KEGG analysis
ekegg <- clusterProfiler::enrichKEGG(gene = entrez,
                                     organism = 'mmu',
                                     pvalueCutoff = 0.05)
svglite(paste0(project_name, "_DB/", project_name,
               "_H89_true_peaks_clusterProfiler_KEGG_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(ekegg)
dev.off()

# Compare KEGG analysis
compKEGG <- clusterProfiler::compareCluster(geneCluster = genes, 
                           fun = "enrichKEGG",
                           organism = "mouse",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_true_peaks_clusterProfiler_KEGG_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(compKEGG, showCategory = 20,
        title = "KEGG Pathway Enrichment Analysis") +
    ggplot2::scale_y_discrete(labels=function(x) str_wrap(x, width=70))
dev.off()
```

#### 11b. P300

```R
project_name <- "P300"

loadGRanges <- function(file_name) {
    tmp <- fread(file_name)
    colnames(tmp) <- c("chr", "start", "end", "name", "score", "strand")
    tmp2 <- makeGRangesFromDataFrame(tmp, keep.extra.columns=TRUE)
    seqlevelsStyle(tmp2) <- "UCSC"
    return(tmp2)
}

# Load data
samplefiles <- list.files(paste0(project_name, "_DB/"),
                          pattern= "_true_peaks.bed", full.names=T)
samplefiles <- sapply(file.path(samplefiles), loadGRanges)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("DMSO", "H89")

txdb         <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)

fwrite(as.data.table(peakAnnoList$DMSO@annoStat),
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_H89_feature_distribution.csv"))

fwrite(as.data.table(peakAnnoList$H89@annoStat),
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_DMSO_feature_distribution.csv"))

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_H89-v-DMSO_feature_distribution_",
               Sys.Date(), ".svg"))
plotAnnoBar(peakAnnoList)
dev.off()

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_H89-v-DMSO_TSS_distribution_",
               Sys.Date(), ".svg"))
plotDistToTSS(peakAnnoList,
    title="Distribution of transcription factor-binding loci \n relative to TSS")
dev.off()

# H89 treated samples
# As a GRanges object, just the peak annotations for H89 for example
P300_H89_annot <- data.frame(peakAnnoList[["H89"]]@anno)

# Get the entrez IDs
entrez_IDs <- P300_H89_annot$geneId

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                                         keys = entrez_IDs,
                                         columns = c("GENENAME"),
                                         keytype = "ENTREZID")

# Change IDs to character type to merge
annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)

# Write to file
P300_H89_annot %>% 
  left_join(annotations_edb, by=c("geneId"="ENTREZID")) %>% 
  write.table(file=paste0(project_name, "_DB/", project_name,
                          "_H89_peak_annotation.txt"),
              sep="\t", quote=F, row.names=F)
  
# Run clusterProfiler's GO enrichment analysis 
ego <- clusterProfiler::enrichGO(gene = entrez_IDs, 
                                 keyType = "ENTREZID", 
                                 OrgDb = org.Mm.eg.db, 
                                 ont = "BP", 
                                 pAdjustMethod = "BH", 
                                 qvalueCutoff = 0.05, 
                                 readable = TRUE)
# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary,
          paste0(project_name, "_DB/", project_name,
                 "_H89_true_peaks_clusterProfiler_GO.csv"))

# Dotplot visualization
svglite(paste0(project_name, "_DB/", project_name,
               "_H89_true_peaks_clusterProfiler_GO_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(ego, showCategory=50) +
    ggplot2::scale_y_discrete(labels=function(x) str_wrap(x, width=70))
dev.off()

# Compare GO analysis
## Create a list with genes from each sample
genes  <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
compGO <- clusterProfiler::compareCluster(geneCluster = genes, 
                                          fun = "enrichGO",
                                          OrgDb = org.Mm.eg.db, 
                                          pvalueCutoff  = 0.05, 
                                          pAdjustMethod = "BH")
svglite(paste0(project_name, "_DB/", project_name,
                "_H89-v-DMSO_true_peaks_clusterProfiler_GO_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(compGO, showCategory = 20,
        title = "GO Pathway Enrichment Analysis") +
    ggplot2::scale_y_discrete(labels=function(x) str_wrap(x, width=70))
dev.off()

# Evaluate and compare KEGG analysis
ekegg <- clusterProfiler::enrichKEGG(gene = entrez,
                                     organism = 'mmu',
                                     pvalueCutoff = 0.05)
svglite(paste0(project_name, "_DB/", project_name,
               "_H89_true_peaks_clusterProfiler_KEGG_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(ekegg)
dev.off()

# Compare KEGG analysis
compKEGG <- clusterProfiler::compareCluster(geneCluster = genes, 
                           fun = "enrichKEGG",
                           organism = "mouse",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_true_peaks_clusterProfiler_KEGG_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(compKEGG, showCategory = 20,
        title = "KEGG Pathway Enrichment Analysis") +
    ggplot2::scale_y_discrete(labels=function(x) str_wrap(x, width=70))
dev.off()

# DMSO treated samples
# As a GRanges object, just the peak annotations for DMSO for example
P300_DMSO_annot <- data.frame(peakAnnoList[["DMSO"]]@anno)

# Get the entrez IDs
entrez_IDs <- P300_DMSO_annot$geneId

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                                         keys = entrez_IDs,
                                         columns = c("GENENAME"),
                                         keytype = "ENTREZID")

# Change IDs to character type to merge
annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)

# Write to file
P300_DMSO_annot %>% 
  left_join(annotations_edb, by=c("geneId"="ENTREZID")) %>% 
  write.table(file=paste0(project_name, "_DB/", project_name,
                          "_DMSO_peak_annotation.txt"),
              sep="\t", quote=F, row.names=F)
  
# Run clusterProfiler's GO enrichment analysis 
ego <- clusterProfiler::enrichGO(gene = entrez_IDs, 
                                 keyType = "ENTREZID", 
                                 OrgDb = org.Mm.eg.db, 
                                 ont = "BP", 
                                 pAdjustMethod = "BH", 
                                 qvalueCutoff = 0.05, 
                                 readable = TRUE)
# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary,
          paste0(project_name, "_DB/", project_name,
                 "_DMSO_true_peaks_clusterProfiler_GO.csv"))

# Dotplot visualization
svglite(paste0(project_name, "_DB/", project_name,
               "_DMSO_true_peaks_clusterProfiler_GO_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(ego, showCategory=50) +
    ggplot2::scale_y_discrete(labels=function(x) str_wrap(x, width=70))
dev.off()

# Compare GO analysis
## Create a list with genes from each sample
genes  <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
compGO <- clusterProfiler::compareCluster(geneCluster = genes, 
                                          fun = "enrichGO",
                                          OrgDb = org.Mm.eg.db, 
                                          pvalueCutoff  = 0.05, 
                                          pAdjustMethod = "BH")
svglite(paste0(project_name, "_DB/", project_name,
                "_H89-v-DMSO_true_peaks_clusterProfiler_GO_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(compGO, showCategory = 20,
        title = "GO Pathway Enrichment Analysis") +
    ggplot2::scale_y_discrete(labels=function(x) str_wrap(x, width=70))
dev.off()

# Evaluate and compare KEGG analysis
ekegg <- clusterProfiler::enrichKEGG(gene = entrez,
                                     organism = 'mmu',
                                     pvalueCutoff = 0.05)
svglite(paste0(project_name, "_DB/", project_name,
               "_DMSO_true_peaks_clusterProfiler_KEGG_",
               Sys.Date(), ".svg"))
clusterProfiler::dotplot(ekegg)
dev.off()
```