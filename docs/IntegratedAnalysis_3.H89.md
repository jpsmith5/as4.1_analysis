# Integrated Analyses

Isolate the enriched P300 and H3K27Ac regions from treated and control

Overlap those with the ATAC-seq dynamic peaks
 - overlap those with the ATAC-seq cluster1 (down) regions
 - overlap those with the ATAC-seq cluster2 (up) regions
 
What are the targets of these regulatory regions?

Overall goal: 
 1. Which genes from the scRNA seq are changed?
 
  - Have CTRL and TRT marker genes for each scRNA-seq cluster (and for the entire dataset)
 
 2. Which regions from the ATACseq are changed?
 
  - Have the dynamic peaks, and the dynamic peaks that cluster to either INCREASE or DECREASE in accessibility following H89 treatment.
 
 3. How does P300 binding change?
 4. How does H3K27Ac change?
 
 **Link each of these individual components together to identify sets of pathways that change following H89 treatment and the corresponding reduction of renin expression**

## 2. Integrate  H89-treated scRNA-seq with ATAC-seq

Combine the H89 results with the scRNA-seq results to identify specific subsets of genomic regions and genes most affected by H89 treatment. 

### 3a. Copy files for H89 comparisons

1. Copy over RNAseq files
```console
cd $PROCESSED/H89/integration_analysis/

cp $PROCESSED/scRNA/R_analysis/As4.1_H89_TRT_markers.csv .
cp $PROCESSED/scRNA/R_analysis/As4.1_H89_CTRL_markers.csv .
cp $PROCESSED/scRNA/R_analysis/markers/*.csv .
```

2. Copy over ATACseq files
```console
export ATACDIR="$PROCESSED/H89/bam_files/FIMO/tomtom/fimo_composites/"
cp cluster2_H89down.bed6 ./dynamic_peaks_DECREASED.bed
cp cluster1_H89up.bed6 ./dynamic_peaks_INCREASED.bed
cp dynamic_peaks.bed ./dynamic_peaks_ALL.bed
cp ${ATACDIR}/nondynamic_peaks.bed ./nondymamic_peaks.bed
cp ${ATACDIR}/all_peaks.bed ./all_peaks.bed
```

### 3b. Load data

```
library(dplyr)
library(patchwork)
library(svglite)
library(ggplot2)
library(sctransform)
library(Seurat)
library(glmGamPoi)
library(data.table)
library(harmony)
library(fastSave)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(GenomicRanges)
library(fastSave)

directory  <- "$PROCESSED/H89/bam_files/"
tomtom_dir <- "$PROCESSED/H89/bam_files/FIMO/tomtom/"
integration_dir <- "$PROCESSED/H89/integration_analysis/"

custom_theme <- function(base_family = "sans", ...){
  theme_classic(base_family = base_family, base_size = 8, ...) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    aspect.ratio = 1,
    legend.position = "right",
    legend.justification="center",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )
}

set.seed(99)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
scrna_dir <- "$PROCESSED/scRNA/R_analysis/"
```

### 3c. Annotate 

```R
trt_markers  <- fread('As4.1_H89_TRT_markers.csv')
ctrl_markers <- fread('As4.1_H89_CTRL_markers.csv')
all_pos_trt_markers  <- fread('As4.1_H89_all+_h89_markers.csv')
all_pos_ctrl_markers <- fread('As4.1_H89_all+_ctrl_markers.csv')
sec24a_arid3b_pos_trt_markers  <- fread('As4.1_H89_Sec24a-Arid3b+_h89_markers.csv')
sec24a_arid3b_pos_ctrl_markers <- fread('As4.1_H89_Sec24a-Arid3b+_ctrl_markers.csv')
ren1_lo_trt_markers  <- fread('As4.1_H89_Ren1_LO_h89_markers.csv')
ren1_lo_ctrl_markers <- fread('As4.1_H89_Ren1_LO_ctrl_markers.csv')
khdrbs3_acaa1b_pos_trt_markers  <- fread('As4.1_H89_Khdrbs3-Acaa1b+_h89_markers.csv')
khdrbs3_acaa1b_pos_ctrl_markers <- fread('As4.1_H89_Khdrbs3-Acaa1b+_ctrl_markers.csv')
remodel_trt_markers   <- fread('As4.1_H89_remodel_h89_markers.csv')
remodel_ctrl_markers  <- fread('As4.1_H89_remodel_ctrl_markers.csv')
gbp7_pos_trt_markers  <- fread('As4.1_H89_Gbp7+_h89_markers.csv')
gbp7_pos_ctrl_markers <- fread('As4.1_H89_Gbp7+_ctrl_markers.csv')

listColumns(EnsDb.Mmusculus.v79)

marker_files <- list("trt_markers", "ctrl_markers",
                     "all_pos_trt_markers", "all_pos_ctrl_markers",
                     "sec24a_arid3b_pos_trt_markers", "sec24a_arid3b_pos_ctrl_markers",
                     "ren1_lo_trt_markers", "ren1_lo_ctrl_markers",
                     "khdrbs3_acaa1b_pos_trt_markers", "khdrbs3_acaa1b_pos_ctrl_markers",
                     "remodel_trt_markers", "remodel_ctrl_markers",
                     "gbp7_pos_trt_markers", "gbp7_pos_ctrl_markers")
annotateMarkers <- function(dat) {
    dt <- get(dat)
    gene_ids <- ensembldb::select(EnsDb.Mmusculus.v79,
        keys= dt[["gene"]],
        keytype = "SYMBOL",
        columns = c("SYMBOL","GENEID", "SEQNAME", "GENESEQSTART",
                    "GENESEQEND", "SEQSTRAND"))
    by  <- join_by(gene == SYMBOL)
    return(dplyr::inner_join(dt, gene_ids, by))
}
                     
marker_files_annotated <- lapply(marker_files, annotateMarkers)
names(marker_files_annotated) <- c("trt_markers", "ctrl_markers",
                     "all_pos_trt_markers", "all_pos_ctrl_markers",
                     "sec24a_arid3b_pos_trt_markers", "sec24a_arid3b_pos_ctrl_markers",
                     "ren1_lo_trt_markers", "ren1_lo_ctrl_markers",
                     "khdrbs3_acaa1b_pos_trt_markers", "khdrbs3_acaa1b_pos_ctrl_markers",
                     "remodel_trt_markers", "remodel_ctrl_markers",
                     "gbp7_pos_trt_markers", "gbp7_pos_ctrl_markers")

for (file in names(marker_files_annotated)) {
    fwrite(marker_files_annotated[[file]], paste0("As4.1_H89_annotated_", file, ".csv"), sep=",")
}

marker_files_annotated_GR <- lapply(marker_files_annotated,
    makeGRangesFromDataFrame, keep.extra.columns=T)
marker_files_annotated_GR <-  GRangesList(marker_files_annotated_GR)
seqlevelsStyle(marker_files_annotated_GR) <- "UCSC"
```

### 3d. Link ATAC-seq to scRNA-seq data 

Take gene expression data.  Link the ATAC-seq data to the nearest expressed gene.

```R
dynamic_decreased <- fread('dynamic_peaks_DECREASED.bed',
                            col.names = c("chr", "start", "end",
                                          "name", "score", "strand"))
dynamic_increased <- fread('dynamic_peaks_INCREASED.bed',
                            col.names = c("chr", "start", "end",
                                          "name", "score", "strand"))
dynamic_all <- fread('dynamic_peaks_ALL.bed',
                      col.names = c("chr", "start", "end"))
nondynamic  <- fread('nondymamic_peaks.bed',
                      col.names = c("chr", "start", "end"))
all_peaks   <- fread('all_peaks.bed',
                     col.names = c("chr", "start", "end"))

dyn_dec_GR <- makeGRangesFromDataFrame(dynamic_decreased, keep.extra.columns=T)
dyn_inc_GR <- makeGRangesFromDataFrame(dynamic_increased, keep.extra.columns=T)
dyn_all_GR <- makeGRangesFromDataFrame(dynamic_all, keep.extra.columns=T)
nondyn_GR  <- makeGRangesFromDataFrame(nondynamic, keep.extra.columns=T)
all_peaks_GR <- makeGRangesFromDataFrame(all_peaks, keep.extra.columns=T)

dtn_inc <- distanceToNearest(x = dyn_inc_GR,
    subject = marker_files_annotated_GR[["ctrl_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
inc_REs <- as.data.table(base::merge(as.data.frame(dyn_inc_GR), 
    as.data.frame(marker_files_annotated_GR[["ctrl_markers"]][subjectHits(dtn_inc),]),
    by.x = 0, by.y = 0))
    
colnames(inc_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_inc_GR`)
## Split `name` to get sort order
inc_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
inc_REs$cluster_no <- as.numeric(inc_REs$cluster_no)
inc_REs <- inc_REs[order(cluster_no),]
## Now add in the distance since this order matches the hit object (`dtn_inc`)
inc_REs[ ,distanceToNearest := mcols(dtn_inc)$distance]
fwrite(inc_REs, "H89_increased-accessibilty-regions_linked_scRNAseq.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(inc_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "H89_increased-accessibilty-regions_linked_scRNAseq_table.csv", sep=",")
```

### 3e. Perform GSEA on these regions/genes

```R
GO_dir <- "$PROCESSED/H89/integration_analysis/GSEA_MSigDB/"
GO_file_hallmarks          <- file.path(GO_dir, "mh.all.v2023.2.Mm.symbols.gmt")
GO_file_canonical_pathways <- file.path(GO_dir, "m2.cp.v2023.2.Mm.symbols.gmt")
GO_file_cell_sig           <- file.path(GO_dir, "m8.all.v2023.2.Mm.symbols.gmt")
GO_file_CGP                <- file.path(GO_dir, "m2.cgp.v2023.2.Mm.symbols.gmt")
GO_file_GTRD  <- file.path(GO_dir, "m3.gtrd.v2023.2.Mm.symbols.gmt")
GO_file_miRNA <- file.path(GO_dir, "m3.mirdb.v2023.2.Mm.symbols.gmt")
GO_file_BP    <- file.path(GO_dir, "m5.go.bp.v2023.2.Mm.symbols.gmt")
GO_file_CC    <- file.path(GO_dir, "m5.go.cc.v2023.2.Mm.symbols.gmt")
GO_file_MF    <- file.path(GO_dir, "m5.go.mf.v2023.2.Mm.symbols.gmt")

GO_list <- c(GO_file_hallmarks, GO_file_canonical_pathways, 
             GO_file_cell_sig, GO_file_CGP, GO_file_GTRD, GO_file_miRNA,
             GO_file_BP, GO_file_CC, GO_file_MF)
names(GO_list) <- c("hallmarks", "canonical", "cell-signature",
                    "CGP", "GTRD", "miRNA", "BP", "CC", "MF")
```

 - Perform GSEA on the dynamic *increased* ATAC regions that overlap RNA-seq markers
```R
inc_REs
inc_REs_gene_list = inc_REs$avg_log2FC
names(inc_REs_gene_list) = inc_REs$gene_symbol
inc_REs_gene_list = sort(inc_REs_gene_list, decreasing = TRUE)
inc_REs_gene_list = inc_REs_gene_list[!duplicated(names(inc_REs_gene_list))]
head(inc_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(inc_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("H89_increased-accessibilty-regions_linked_scRNAseq_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}
```

### 3f. RNA marker list versus decreased dynamic ATAC-seq peaks

```R
dtn_dec <- distanceToNearest(x = dyn_dec_GR,
    subject = marker_files_annotated_GR[["ctrl_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
dec_REs <- as.data.table(base::merge(as.data.frame(dyn_dec_GR), 
    as.data.frame(marker_files_annotated_GR[["ctrl_markers"]][subjectHits(dtn_dec),]),
    by.x = 0, by.y = 0))
    
colnames(dec_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_dec_GR`)
## Split `name` to get sort order
dec_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
dec_REs$cluster_no <- as.numeric(dec_REs$cluster_no)
dec_REs <- dec_REs[order(cluster_no),]
## Now add in the distance sdece this order matches the hit object (`dtn_dec`)
dec_REs[ ,distanceToNearest := mcols(dtn_dec)$distance]
fwrite(dec_REs, "H89_decreased-accessibilty-regions_linked_scRNAseq.csv", sep=",")
```

Which genes have the most REs?

```R
gene_RE_table <- as.data.table(sort(table(dec_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "H89_decreased-accessibilty-regions_linked_scRNAseq_table.csv", sep=",")
```

 - Perform GSEA on the dynamic *decreased* ATAC regions that overlap RNA-seq markers
```R
dec_REs
dec_REs_gene_list = dec_REs$avg_log2FC
names(dec_REs_gene_list) = dec_REs$gene_symbol
dec_REs_gene_list = sort(dec_REs_gene_list, decreasing = TRUE)
dec_REs_gene_list = dec_REs_gene_list[!duplicated(names(dec_REs_gene_list))]
head(dec_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(dec_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("H89_decreased-accessibilty-regions_linked_scRNAseq_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}
```

This overlapping gene set between increased and decreased regions is:
```R
ovp_genes <- as.data.table(sort(dec_REs_gene_list[dec_REs_gene_list %in% inc_REs_gene_list]))
colnames(ovp_genes) <- "avg_log2FC"
ovp_genes[, gene_symbol := names(sort(dec_REs_gene_list[dec_REs_gene_list %in% inc_REs_gene_list]))]
fwrite(ovp_genes, "H89_decreased-accessibilty-regions_linked_scRNAseq_IN_increased-accessibility-gene-list.csv", sep=",")
```

### 3g. Evaluate open chromatin regions and gene expression overlaps in each scRNA-seq population

1. all+
2. remodel
3. Khdrbs3|Acaa1b+
4. Ren1_LO
5. Sec24a|Arid3b+
6. Gbp7+

#### "all+" population

```R
dtn_inc <- distanceToNearest(x = dyn_inc_GR,
    subject = marker_files_annotated_GR[["all_pos_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
inc_REs <- as.data.table(base::merge(as.data.frame(dyn_inc_GR), 
    as.data.frame(marker_files_annotated_GR[["all_pos_trt_markers"]][subjectHits(dtn_inc),]),
    by.x = 0, by.y = 0))
    
colnames(inc_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_inc_GR`)
## Split `name` to get sort order
inc_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
inc_REs$cluster_no <- as.numeric(inc_REs$cluster_no)
inc_REs <- inc_REs[order(cluster_no),]
## Now add in the distance since this order matches the hit object (`dtn_inc`)
inc_REs[ ,distanceToNearest := mcols(dtn_inc)$distance]
fwrite(inc_REs, "dynamic_increased_all-pos_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(inc_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_increased_all-pos_linked-REs_table.csv", sep=",")

inc_REs
inc_REs_gene_list = inc_REs$avg_log2FC
names(inc_REs_gene_list) = inc_REs$gene_symbol
inc_REs_gene_list = sort(inc_REs_gene_list, decreasing = TRUE)
inc_REs_gene_list = inc_REs_gene_list[!duplicated(names(inc_REs_gene_list))]
head(inc_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(inc_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_increased_all-pos_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

dtn_dec <- distanceToNearest(x = dyn_dec_GR,
    subject = marker_files_annotated_GR[["all_pos_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
dec_REs <- as.data.table(base::merge(as.data.frame(dyn_dec_GR), 
    as.data.frame(marker_files_annotated_GR[["all_pos_trt_markers"]][subjectHits(dtn_dec),]),
    by.x = 0, by.y = 0))
    
colnames(dec_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_dec_GR`)
## Split `name` to get sort order
dec_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
dec_REs$cluster_no <- as.numeric(dec_REs$cluster_no)
dec_REs <- dec_REs[order(cluster_no),]
## Now add in the distance sdece this order matches the hit object (`dtn_dec`)
dec_REs[ ,distanceToNearest := mcols(dtn_dec)$distance]
fwrite(dec_REs, "dynamic_decreased_all-pos_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(dec_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_decreased_all-pos_linked-REs_table.csv", sep=",")

dec_REs
dec_REs_gene_list = dec_REs$avg_log2FC
names(dec_REs_gene_list) = dec_REs$gene_symbol
dec_REs_gene_list = sort(dec_REs_gene_list, decreasing = TRUE)
dec_REs_gene_list = dec_REs_gene_list[!duplicated(names(dec_REs_gene_list))]
head(dec_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(dec_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_decreased_all-pos_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}
```

#### "remodel" population

```R
dtn_inc <- distanceToNearest(x = dyn_inc_GR,
    subject = marker_files_annotated_GR[["remodel_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
inc_REs <- as.data.table(base::merge(as.data.frame(dyn_inc_GR), 
    as.data.frame(marker_files_annotated_GR[["remodel_trt_markers"]][subjectHits(dtn_inc),]),
    by.x = 0, by.y = 0))
    
colnames(inc_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_inc_GR`)
## Split `name` to get sort order
inc_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
inc_REs$cluster_no <- as.numeric(inc_REs$cluster_no)
inc_REs <- inc_REs[order(cluster_no),]
## Now add in the distance since this order matches the hit object (`dtn_inc`)
inc_REs[ ,distanceToNearest := mcols(dtn_inc)$distance]
fwrite(inc_REs, "dynamic_increased_remodel_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(inc_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_increased_remodel_linked-REs_table.csv", sep=",")

inc_REs
inc_REs_gene_list = inc_REs$avg_log2FC
names(inc_REs_gene_list) = inc_REs$gene_symbol
inc_REs_gene_list = sort(inc_REs_gene_list, decreasing = TRUE)
inc_REs_gene_list = inc_REs_gene_list[!duplicated(names(inc_REs_gene_list))]
head(inc_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(inc_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_increased_remodel_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

dtn_dec <- distanceToNearest(x = dyn_dec_GR,
    subject = marker_files_annotated_GR[["remodel_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
dec_REs <- as.data.table(base::merge(as.data.frame(dyn_dec_GR), 
    as.data.frame(marker_files_annotated_GR[["remodel_trt_markers"]][subjectHits(dtn_dec),]),
    by.x = 0, by.y = 0))
    
colnames(dec_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_dec_GR`)
## Split `name` to get sort order
dec_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
dec_REs$cluster_no <- as.numeric(dec_REs$cluster_no)
dec_REs <- dec_REs[order(cluster_no),]
## Now add in the distance sdece this order matches the hit object (`dtn_dec`)
dec_REs[ ,distanceToNearest := mcols(dtn_dec)$distance]
fwrite(dec_REs, "dynamic_decreased_remodel_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(dec_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_decreased_remodel_linked-REs_table.csv", sep=",")

dec_REs
dec_REs_gene_list = dec_REs$avg_log2FC
names(dec_REs_gene_list) = dec_REs$gene_symbol
dec_REs_gene_list = sort(dec_REs_gene_list, decreasing = TRUE)
dec_REs_gene_list = dec_REs_gene_list[!duplicated(names(dec_REs_gene_list))]
head(dec_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(dec_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_decreased_remodel_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}
```

#### "Khdrbs3|Acaa1b+" population

```R
dtn_inc <- distanceToNearest(x = dyn_inc_GR,
    subject = marker_files_annotated_GR[["khdrbs3_acaa1b_pos_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
inc_REs <- as.data.table(base::merge(as.data.frame(dyn_inc_GR), 
    as.data.frame(marker_files_annotated_GR[["khdrbs3_acaa1b_pos_trt_markers"]][subjectHits(dtn_inc),]),
    by.x = 0, by.y = 0))
    
colnames(inc_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_inc_GR`)
## Split `name` to get sort order
inc_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
inc_REs$cluster_no <- as.numeric(inc_REs$cluster_no)
inc_REs <- inc_REs[order(cluster_no),]
## Now add in the distance since this order matches the hit object (`dtn_inc`)
inc_REs[ ,distanceToNearest := mcols(dtn_inc)$distance]
fwrite(inc_REs, "dynamic_increased_khdrbs3_acaa1b_pos_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(inc_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_increased_khdrbs3_acaa1b_pos_linked-REs_table.csv", sep=",")

inc_REs
inc_REs_gene_list = inc_REs$avg_log2FC
names(inc_REs_gene_list) = inc_REs$gene_symbol
inc_REs_gene_list = sort(inc_REs_gene_list, decreasing = TRUE)
inc_REs_gene_list = inc_REs_gene_list[!duplicated(names(inc_REs_gene_list))]
head(inc_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(inc_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_increased_khdrbs3_acaa1b_pos_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

dtn_dec <- distanceToNearest(x = dyn_dec_GR,
    subject = marker_files_annotated_GR[["khdrbs3_acaa1b_pos_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
dec_REs <- as.data.table(base::merge(as.data.frame(dyn_dec_GR), 
    as.data.frame(marker_files_annotated_GR[["khdrbs3_acaa1b_pos_trt_markers"]][subjectHits(dtn_dec),]),
    by.x = 0, by.y = 0))
    
colnames(dec_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_dec_GR`)
## Split `name` to get sort order
dec_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
dec_REs$cluster_no <- as.numeric(dec_REs$cluster_no)
dec_REs <- dec_REs[order(cluster_no),]
## Now add in the distance sdece this order matches the hit object (`dtn_dec`)
dec_REs[ ,distanceToNearest := mcols(dtn_dec)$distance]
fwrite(dec_REs, "dynamic_decreased_khdrbs3_acaa1b_pos_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(dec_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_decreased_khdrbs3_acaa1b_pos_linked-REs_table.csv", sep=",")

dec_REs
dec_REs_gene_list = dec_REs$avg_log2FC
names(dec_REs_gene_list) = dec_REs$gene_symbol
dec_REs_gene_list = sort(dec_REs_gene_list, decreasing = TRUE)
dec_REs_gene_list = dec_REs_gene_list[!duplicated(names(dec_REs_gene_list))]
head(dec_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(dec_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_decreased_khdrbs3_acaa1b_pos_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}
```

#### "Ren1_LO" population

```R
dtn_inc <- distanceToNearest(x = dyn_inc_GR,
    subject = marker_files_annotated_GR[["ren1_lo_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
inc_REs <- as.data.table(base::merge(as.data.frame(dyn_inc_GR), 
    as.data.frame(marker_files_annotated_GR[["ren1_lo_trt_markers"]][subjectHits(dtn_inc),]),
    by.x = 0, by.y = 0))
    
colnames(inc_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_inc_GR`)
## Split `name` to get sort order
inc_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
inc_REs$cluster_no <- as.numeric(inc_REs$cluster_no)
inc_REs <- inc_REs[order(cluster_no),]
## Now add in the distance since this order matches the hit object (`dtn_inc`)
inc_REs[ ,distanceToNearest := mcols(dtn_inc)$distance]
fwrite(inc_REs, "dynamic_increased_ren1-lo_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(inc_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_increased_ren1-lo_linked-REs_table.csv", sep=",")

inc_REs
inc_REs_gene_list = inc_REs$avg_log2FC
names(inc_REs_gene_list) = inc_REs$gene_symbol
inc_REs_gene_list = sort(inc_REs_gene_list, decreasing = TRUE)
inc_REs_gene_list = inc_REs_gene_list[!duplicated(names(inc_REs_gene_list))]
head(inc_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(inc_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_increased_ren1-lo_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

dtn_dec <- distanceToNearest(x = dyn_dec_GR,
    subject = marker_files_annotated_GR[["ren1_lo_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
dec_REs <- as.data.table(base::merge(as.data.frame(dyn_dec_GR), 
    as.data.frame(marker_files_annotated_GR[["ren1_lo_trt_markers"]][subjectHits(dtn_dec),]),
    by.x = 0, by.y = 0))
    
colnames(dec_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_dec_GR`)
## Split `name` to get sort order
dec_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
dec_REs$cluster_no <- as.numeric(dec_REs$cluster_no)
dec_REs <- dec_REs[order(cluster_no),]
## Now add in the distance sdece this order matches the hit object (`dtn_dec`)
dec_REs[ ,distanceToNearest := mcols(dtn_dec)$distance]
fwrite(dec_REs, "dynamic_decreased_ren1-lo_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(dec_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_decreased_ren1-lo_linked-REs_table.csv", sep=",")

dec_REs
dec_REs_gene_list = dec_REs$avg_log2FC
names(dec_REs_gene_list) = dec_REs$gene_symbol
dec_REs_gene_list = sort(dec_REs_gene_list, decreasing = TRUE)
dec_REs_gene_list = dec_REs_gene_list[!duplicated(names(dec_REs_gene_list))]
head(dec_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(dec_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_decreased_ren1-lo_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}
```

#### "Sec24a|Arid3b+" population

```R
dtn_inc <- distanceToNearest(x = dyn_inc_GR,
    subject = marker_files_annotated_GR[["sec24a_arid3b_pos_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
inc_REs <- as.data.table(base::merge(as.data.frame(dyn_inc_GR), 
    as.data.frame(marker_files_annotated_GR[["sec24a_arid3b_pos_trt_markers"]][subjectHits(dtn_inc),]),
    by.x = 0, by.y = 0))
    
colnames(inc_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_inc_GR`)
## Split `name` to get sort order
inc_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
inc_REs$cluster_no <- as.numeric(inc_REs$cluster_no)
inc_REs <- inc_REs[order(cluster_no),]
## Now add in the distance since this order matches the hit object (`dtn_inc`)
inc_REs[ ,distanceToNearest := mcols(dtn_inc)$distance]
fwrite(inc_REs, "dynamic_increased_sec24a-arid3b-pos_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(inc_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_increased_sec24a-arid3b-pos_linked-REs_table.csv", sep=",")

inc_REs
inc_REs_gene_list = inc_REs$avg_log2FC
names(inc_REs_gene_list) = inc_REs$gene_symbol
inc_REs_gene_list = sort(inc_REs_gene_list, decreasing = TRUE)
inc_REs_gene_list = inc_REs_gene_list[!duplicated(names(inc_REs_gene_list))]
head(inc_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(inc_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_increased_sec24a-arid3b-pos_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

dtn_dec <- distanceToNearest(x = dyn_dec_GR,
    subject = marker_files_annotated_GR[["sec24a_arid3b_pos_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
dec_REs <- as.data.table(base::merge(as.data.frame(dyn_dec_GR), 
    as.data.frame(marker_files_annotated_GR[["sec24a_arid3b_pos_trt_markers"]][subjectHits(dtn_dec),]),
    by.x = 0, by.y = 0))
    
colnames(dec_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_dec_GR`)
## Split `name` to get sort order
dec_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
dec_REs$cluster_no <- as.numeric(dec_REs$cluster_no)
dec_REs <- dec_REs[order(cluster_no),]
## Now add in the distance sdece this order matches the hit object (`dtn_dec`)
dec_REs[ ,distanceToNearest := mcols(dtn_dec)$distance]
fwrite(dec_REs, "dynamic_decreased_sec24a-arid3b-pos_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(dec_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_decreased_sec24a-arid3b-pos_linked-REs_table.csv", sep=",")

dec_REs
dec_REs_gene_list = dec_REs$avg_log2FC
names(dec_REs_gene_list) = dec_REs$gene_symbol
dec_REs_gene_list = sort(dec_REs_gene_list, decreasing = TRUE)
dec_REs_gene_list = dec_REs_gene_list[!duplicated(names(dec_REs_gene_list))]
head(dec_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(dec_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_decreased_sec24a-arid3b-pos_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}
```

#### "Gbp7+" population

```R
dtn_inc <- distanceToNearest(x = dyn_inc_GR,
    subject = marker_files_annotated_GR[["gbp7_pos_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
inc_REs <- as.data.table(base::merge(as.data.frame(dyn_inc_GR), 
    as.data.frame(marker_files_annotated_GR[["gbp7_pos_trt_markers"]][subjectHits(dtn_inc),]),
    by.x = 0, by.y = 0))
    
colnames(inc_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_inc_GR`)
## Split `name` to get sort order
inc_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
inc_REs$cluster_no <- as.numeric(inc_REs$cluster_no)
inc_REs <- inc_REs[order(cluster_no),]
## Now add in the distance since this order matches the hit object (`dtn_inc`)
inc_REs[ ,distanceToNearest := mcols(dtn_inc)$distance]
fwrite(inc_REs, "dynamic_increased_gbp7-pos_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(inc_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_increased_gbp7-pos_linked-REs_table.csv", sep=",")

inc_REs
inc_REs_gene_list = inc_REs$avg_log2FC
names(inc_REs_gene_list) = inc_REs$gene_symbol
inc_REs_gene_list = sort(inc_REs_gene_list, decreasing = TRUE)
inc_REs_gene_list = inc_REs_gene_list[!duplicated(names(inc_REs_gene_list))]
head(inc_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(inc_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_increased_gbp7-pos_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

dtn_dec <- distanceToNearest(x = dyn_dec_GR,
    subject = marker_files_annotated_GR[["gbp7_pos_trt_markers"]],
    ignore.strand=FALSE)

# Merge gene expression with atac-seq
dec_REs <- as.data.table(base::merge(as.data.frame(dyn_dec_GR), 
    as.data.frame(marker_files_annotated_GR[["gbp7_pos_trt_markers"]][subjectHits(dtn_dec),]),
    by.x = 0, by.y = 0))
    
colnames(dec_REs) <- c("row_names", "chr_atac", "start_atac", "end_atac", "width_atac", "strand_atac", "name", "score", "chr_rna", "start_rna", "end_rna", "width_rna", "strand_rna", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene_symbol", "gene_id", "seq_strand")

# Add actual distance to features
## First, must sort the object to match the ATAC-seq object (`dyn_dec_GR`)
## Split `name` to get sort order
dec_REs[, c("cluster_name", "cluster_no") := tstrsplit(name, "_", fixed=TRUE)]
dec_REs$cluster_no <- as.numeric(dec_REs$cluster_no)
dec_REs <- dec_REs[order(cluster_no),]
## Now add in the distance sdece this order matches the hit object (`dtn_dec`)
dec_REs[ ,distanceToNearest := mcols(dtn_dec)$distance]
fwrite(dec_REs, "dynamic_decreased_gbp7-pos_linked-REs.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(dec_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "dynamic_decreased_gbp7-pos_linked-REs_table.csv", sep=",")

dec_REs
dec_REs_gene_list = dec_REs$avg_log2FC
names(dec_REs_gene_list) = dec_REs$gene_symbol
dec_REs_gene_list = sort(dec_REs_gene_list, decreasing = TRUE)
dec_REs_gene_list = dec_REs_gene_list[!duplicated(names(dec_REs_gene_list))]
head(dec_REs_gene_list)

for (GO_file in names(GO_list)) {
    res = GSEA(dec_REs_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("dynamic_decreased_gbp7-pos_linked-REs_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}
```

### 3h. Link ChIP-seq data

```R
library(ATACseqQC)
library(MotifDb)
library(DiffBind)
library(svglite)
library(stringr)
library(profileplyr)
library(ChIPpeakAnno)
library(data.table)
library(org.Mm.eg.db)
library(reactome.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(scales)
library(stringi)
library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(dplyr)
library(rtracklayer)
library(fastSave)
library(DESeq2)
library(lattice)
library(bigWig)
library(genefilter)
library(ggplot2)
library(ggrepel)

set.seed <- 99

txdb   <- TxDb.Mmusculus.UCSC.mm10.knownGene
genome <- BSgenome.Mmusculus.UCSC.mm10

.validateInputs = function(checkList) {
    nms = names(checkList)
    for(i in seq_along(checkList)){
        fail = FALSE
        clss = checkList[[i]]
        x = get(nms[i], envir=parent.frame(1))
        for(cls in clss){
            if (is(x, cls)) fail = append(fail, TRUE)
        }
        if(!any(fail)) 
            stop(paste0(nms[i], " must be a ", paste(clss, collapse=" or "), 
                        ".  Got: ", class(x)))
    }
}

# Add raw counts from the signal tracks overlapping a bed file
getRawCountsInterval <- function(df, bigWig_path, file_prefix = 'M') {
    c1 = list(df=c("data.table", "data.frame"))
    .validateInputs(c1)
	df = df[,c("chr", "start", "end")]
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
	rowname_vals <- paste0(df$chr, ':', df$start, '-', df$end)
	row.names(inten_df) <- rowname_vals
	return(inten_df)
}

# Convert a raw counts interval object to a GRanges object
countsToGR <- function(counts_df) {
    dt <- as.data.table(counts_df, keep.rownames=T)
    vec_names <- colnames(dt[,-1])
    colnames(dt) <- c("region", paste0("counts_",
                                       seq(from=1, to=ncol(dt[,-1]), by=1)))
    dt[, c("chr", "coor") := tstrsplit(region, ":", fixed=TRUE)]
    dt[, c("start", "end") := tstrsplit(coor, "-", fixed=TRUE)]
    final_dt <- data.table(chr=dt$chr,
                           start=dt$start,
                           end=dt$end,
                           name=dt$region,
                           counts_1=dt$counts_1,
                           counts_2=dt$counts_2)
    colnames(final_dt) <- c("chr", "start", "end", "name", vec_names)
    return(makeGRangesFromDataFrame(final_dt, keep.extra.columns=TRUE))
}
```

- Load the ATAC all peak files

```R
atac_dir <- "$PROCESSED/H89/bam_files/"
DEseq_region_sets <- c(paste0(atac_dir, 'notDifferent_As4.1_H89-vs-DMSO_1e-04_FDR.bed'), 
                       paste0(atac_dir, 'increased_As4.1_H89-vs-DMSO_1e-04_FDR.bed'), 
                       paste0(atac_dir, 'decreased_As4.1_H89-vs-DMSO_1e-04_FDR.bed'))
                       
consensus_peaks_file <- fread(paste0(atac_dir, "As4.1_smallpval_rmBlacklist.bed"))
colnames(consensus_peaks_file) <- c("chr", "start", "end")

consensus_peaks_tsv <- read.table(paste0(atac_dir, "As4.1_smallpval_rmBlacklist.bed"))
colnames(consensus_peaks_tsv) <- c("chr", "start", "end")
trt_H89_counts   <- getRawCountsInterval(consensus_peaks_tsv, atac_dir, file_prefix = 'H89')
ctrl_DMSO_counts <- getRawCountsInterval(consensus_peaks_tsv, atac_dir, file_prefix = 'DMSO')

trt_GR  <- countsToGR(trt_H89_counts)
ctrl_GR <- countsToGR(ctrl_DMSO_counts)

trt_ens_GR <- copy(trt_GR)
ctrl_ens_GR <- copy(ctrl_GR)
seqlevelsStyle(trt_ens_GR)  <- "Ensembl"
seqlevelsStyle(ctrl_ens_GR) <- "Ensembl"

consensus_GR <- merge(trt_GR, ctrl_GR)
consensus_ens_GR <- merge(trt_ens_GR, ctrl_ens_GR)

dynamic_decreased <- fread('dynamic_peaks_DECREASED.bed',
                            col.names = c("chr", "start", "end",
                                          "name", "score", "strand"))
dynamic_increased <- fread('dynamic_peaks_INCREASED.bed',
                            col.names = c("chr", "start", "end",
                                          "name", "score", "strand"))
dynamic_all <- fread('dynamic_peaks_ALL.bed',
                      col.names = c("chr", "start", "end"))
nondynamic  <- fread('nondymamic_peaks.bed',
                      col.names = c("chr", "start", "end"))
all_peaks   <- fread('all_peaks.bed',
                     col.names = c("chr", "start", "end"))

dyn_dec_GR <- makeGRangesFromDataFrame(dynamic_decreased, keep.extra.columns=T)
dyn_inc_GR <- makeGRangesFromDataFrame(dynamic_increased, keep.extra.columns=T)
dyn_all_GR <- makeGRangesFromDataFrame(dynamic_all, keep.extra.columns=T)
nondyn_GR  <- makeGRangesFromDataFrame(nondynamic, keep.extra.columns=T)
all_peaks_GR <- makeGRangesFromDataFrame(all_peaks, keep.extra.columns=T)

trt_dyn_inc_counts  <- getRawCountsInterval(dynamic_increased, atac_dir, file_prefix = 'H89')
ctrl_dyn_inc_counts <- getRawCountsInterval(dynamic_increased, atac_dir, file_prefix = 'DMSO')
trt_dyn_inc_GR      <- countsToGR(trt_dyn_inc_counts)
ctrl_dyn_inc_GR     <- countsToGR(ctrl_dyn_inc_counts)
trt_dyn_inc_ens_GR  <- copy(trt_dyn_inc_GR)
ctrl_dyn_inc_ens_GR <- copy(ctrl_dyn_inc_GR)
seqlevelsStyle(trt_dyn_inc_ens_GR)  <- "Ensembl"
seqlevelsStyle(ctrl_dyn_inc_ens_GR) <- "Ensembl"
dyninc_GR     <- merge(trt_dyn_inc_GR, ctrl_dyn_inc_GR)
dyninc_ens_GR <- merge(trt_dyn_inc_ens_GR, ctrl_dyn_inc_ens_GR)


trt_dyn_dec_counts  <- getRawCountsInterval(dynamic_decreased, atac_dir, file_prefix = 'H89')
ctrl_dyn_dec_counts <- getRawCountsInterval(dynamic_decreased, atac_dir, file_prefix = 'DMSO')
trt_dyn_dec_GR      <- countsToGR(trt_dyn_dec_counts)
ctrl_dyn_dec_GR     <- countsToGR(ctrl_dyn_dec_counts)
trt_dyn_dec_ens_GR  <- copy(trt_dyn_dec_GR)
ctrl_dyn_dec_ens_GR <- copy(ctrl_dyn_dec_GR)
seqlevelsStyle(trt_dyn_dec_ens_GR)  <- "Ensembl"
seqlevelsStyle(ctrl_dyn_dec_ens_GR) <- "Ensembl"
dyndec_GR     <- merge(trt_dyn_dec_GR, ctrl_dyn_dec_GR)
dyndec_ens_GR <- merge(trt_dyn_dec_ens_GR, ctrl_dyn_dec_ens_GR)

trt_nondyn_counts  <- getRawCountsInterval(nondynamic, atac_dir, file_prefix = 'H89')
ctrl_nondyn_counts <- getRawCountsInterval(nondynamic, atac_dir, file_prefix = 'DMSO')
trt_nondyn_GR      <- countsToGR(trt_nondyn_counts)
ctrl_nondyn_GR     <- countsToGR(ctrl_nondyn_counts)
trt_nondyn_ens_GR  <- copy(trt_nondyn_GR)
ctrl_nondyn_ens_GR <- copy(ctrl_nondyn_GR)
seqlevelsStyle(trt_nondyn_ens_GR)  <- "Ensembl"
seqlevelsStyle(ctrl_nondyn_ens_GR) <- "Ensembl"
nondyn_GR     <- merge(trt_nondyn_GR, ctrl_nondyn_GR)
nondyn_ens_GR <- merge(trt_nondyn_ens_GR, ctrl_nondyn_ens_GR)
``` 

 - Load ChIP-seq data for H3K27Ac

```R
chip_dir <- "$PROCESSED/chipseq/"
promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)
seqlevelsStyle(promoters) <- "Ensembl"

load.pigz(file.path(chip_dir, 'chip_analysis.RData'))

project_name <- "H3K27Ac"

# DBA object: dba_H3K27Ac
# dba.report: dba_H3K27Ac_DB
# annotated dba.report: dba_H3K27Ac_peaksAnno

H3K27Ac_H89down <- dba_H3K27Ac_peaksAnno[dba_H3K27Ac_peaksAnno$Fold > 0,]
H3K27Ac_H89up   <- dba_H3K27Ac_peaksAnno[dba_H3K27Ac_peaksAnno$Fold < 0,]
```

 - Which promoters overlap H3K27Ac regions?
 
```R
H3K27Ac_H89down_promoter_olp <- subsetByOverlaps(H3K27Ac_H89down, promoters)
H3K27Ac_H89up_promoter_olp   <- subsetByOverlaps(H3K27Ac_H89up, promoters)

fwrite(as.data.table(H3K27Ac_H89down_promoter_olp), "H3K27Ac_H89down_promoter_olp.csv", sep=",")
fwrite(as.data.table(H3K27Ac_H89up_promoter_olp), "H3K27Ac_H89up_promoter_olp.csv", sep=",")
```

 - Which ATAC-seq regions overlap H3K27Ac regions?
 
```R
seqlevelsStyle(H3K27Ac_H89down) <- "UCSC"
seqlevelsStyle(H3K27Ac_H89up)   <- "UCSC"

H3K27Ac_H89down_dyn_inc_olp <- subsetByOverlaps(H3K27Ac_H89down, dyn_inc_GR)
H3K27Ac_H89up_dyn_inc_olp   <- subsetByOverlaps(H3K27Ac_H89up, dyn_inc_GR)
fwrite(as.data.table(H3K27Ac_H89down_dyn_inc_olp), "H3K27Ac_H89down_dyn_inc_olp.csv", sep=",")
fwrite(as.data.table(H3K27Ac_H89up_dyn_inc_olp), "H3K27Ac_H89up_dyn_inc_olp.csv", sep=",")

H3K27Ac_H89down_dyn_dec_olp <- subsetByOverlaps(H3K27Ac_H89down, dyn_dec_GR)
H3K27Ac_H89up_dyn_dec_olp   <- subsetByOverlaps(H3K27Ac_H89up, dyn_dec_GR)
fwrite(as.data.table(H3K27Ac_H89down_dyn_dec_olp), "H3K27Ac_H89down_dyn_dec_olp.csv", sep=",")
fwrite(as.data.table(H3K27Ac_H89up_dyn_dec_olp), "H3K27Ac_H89up_dyn_dec_olp.csv", sep=",")
```

 - Of the ATAC-seq regions that overlap H3K27Ac ChIP regions, which of those associated target genes are present in the marker lists from scRNA-seq?
 - Here, look at which ChIP regions overlap the RNAseq markers
 
```R
H3K27Ac_H89down_trt_markers_olp <- subsetByOverlaps(
    H3K27Ac_H89down, marker_files_annotated_GR[["trt_markers"]])
H3K27Ac_H89up_trt_markers_olp   <- subsetByOverlaps(
    H3K27Ac_H89up, marker_files_annotated_GR[["trt_markers"]])
```

 - Here, look at which ChIP-ATAC overlap regions associated target genes are present in the RNAseq marker list
 
```R
# Increased ATAC regions
markers_olp_H3K27Ac_H89down_and_incATAC <- dplyr::intersect(marker_files_annotated_GR[["trt_markers"]]$gene, H3K27Ac_H89down_dyn_inc_olp$symbol)
H3K27Ac_H89down_incATAC_RNA_olp <- marker_files_annotated_GR[["trt_markers"]][marker_files_annotated_GR[["trt_markers"]]$gene %in% markers_olp_H3K27Ac_H89down_and_incATAC,]
fwrite(as.data.table(H3K27Ac_H89down_incATAC_RNA_olp), "H3K27Ac_H89down-incATAC-RNA_olp.csv", sep=",")

markers_olp_H3K27Ac_H89up_and_incATAC <- dplyr::intersect(marker_files_annotated_GR[["trt_markers"]]$gene, H3K27Ac_H89up_dyn_inc_olp$symbol)
H3K27Ac_H89up_incATAC_RNA_olp <- marker_files_annotated_GR[["trt_markers"]][marker_files_annotated_GR[["trt_markers"]]$gene %in% markers_olp_H3K27Ac_H89up_and_incATAC,]
fwrite(as.data.table(H3K27Ac_H89up_incATAC_RNA_olp), "H3K27Ac_H89up-incATAC-RNA_olp.csv", sep=",")

# Decreased ATAC regions
markers_olp_H3K27Ac_H89down_and_decATAC <- dplyr::intersect(marker_files_annotated_GR[["trt_markers"]]$gene, H3K27Ac_H89down_dyn_dec_olp$symbol)
H3K27Ac_H89down_decATAC_RNA_olp <- marker_files_annotated_GR[["trt_markers"]][marker_files_annotated_GR[["trt_markers"]]$gene %in% markers_olp_H3K27Ac_H89down_and_decATAC,]
fwrite(as.data.table(H3K27Ac_H89down_decATAC_RNA_olp), "H3K27Ac_H89down-decATAC-RNA_olp.csv", sep=",")

markers_olp_H3K27Ac_H89up_and_decATAC <- dplyr::intersect(marker_files_annotated_GR[["trt_markers"]]$gene, H3K27Ac_H89up_dyn_dec_olp$symbol)
H3K27Ac_H89up_decATAC_RNA_olp <- marker_files_annotated_GR[["trt_markers"]][marker_files_annotated_GR[["trt_markers"]]$gene %in% markers_olp_H3K27Ac_H89up_and_decATAC,]
fwrite(as.data.table(H3K27Ac_H89up_decATAC_RNA_olp), "H3K27Ac_H89up-decATAC-RNA_olp.csv", sep=",")
```

 - Here, look at which ChIP-ATAC overlap regions associated target genes are present in the `Ren1_LO` cluster's RNAseq marker list

```R
# Increased ATAC regions
markers_olp_H3K27Ac_H89down_and_incATAC <- dplyr::intersect(marker_files_annotated_GR[["ren1_lo_ren1_lo_trt_markers"]]$gene, H3K27Ac_H89down_dyn_inc_olp$symbol)
H3K27Ac_H89down_incATAC_RNA_olp <- marker_files_annotated_GR[["ren1_lo_trt_markers"]][marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene %in% markers_olp_H3K27Ac_H89down_and_incATAC,]
fwrite(as.data.table(H3K27Ac_H89down_incATAC_RNA_olp), "H3K27Ac_H89down-incATAC-Ren1-LO-RNA_olp.csv", sep=",")

markers_olp_H3K27Ac_H89up_and_incATAC <- dplyr::intersect(marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene, H3K27Ac_H89up_dyn_inc_olp$symbol)
H3K27Ac_H89up_incATAC_RNA_olp <- marker_files_annotated_GR[["ren1_lo_trt_markers"]][marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene %in% markers_olp_H3K27Ac_H89up_and_incATAC,]
fwrite(as.data.table(H3K27Ac_H89up_incATAC_RNA_olp), "H3K27Ac_H89up-incATAC-Ren1-LO-RNA_olp.csv", sep=",")

# Decreased ATAC regions
markers_olp_H3K27Ac_H89down_and_decATAC <- dplyr::intersect(marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene, H3K27Ac_H89down_dyn_dec_olp$symbol)
H3K27Ac_H89down_decATAC_RNA_olp <- marker_files_annotated_GR[["ren1_lo_trt_markers"]][marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene %in% markers_olp_H3K27Ac_H89down_and_decATAC,]
fwrite(as.data.table(H3K27Ac_H89down_decATAC_RNA_olp), "H3K27Ac_H89down-decATAC-Ren1-LO-RNA_olp.csv", sep=",")

markers_olp_H3K27Ac_H89up_and_decATAC <- dplyr::intersect(marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene, H3K27Ac_H89up_dyn_dec_olp$symbol)
H3K27Ac_H89up_decATAC_RNA_olp <- marker_files_annotated_GR[["ren1_lo_trt_markers"]][marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene %in% markers_olp_H3K27Ac_H89up_and_decATAC,]
fwrite(as.data.table(H3K27Ac_H89up_decATAC_RNA_olp), "H3K27Ac_H89up-decATAC-Ren1-LO-RNA_olp.csv", sep=",")
```

 - Load ChIP-seq data for P300
 
```R
project_name <- "P300"
P300_H89down <- dba_P300_peaksAnno[dba_P300_peaksAnno$Fold > 0,]
P300_H89up   <- dba_P300_peaksAnno[dba_P300_peaksAnno$Fold < 0,]
```

 - Which promoters overlap P300 regions?
 
```R
P300_H89down_promoter_olp <- subsetByOverlaps(P300_H89down, promoters)
P300_H89up_promoter_olp   <- subsetByOverlaps(P300_H89up, promoters)

fwrite(as.data.table(P300_H89down_promoter_olp), "P300_H89down_promoter_olp.csv", sep=",")
fwrite(as.data.table(P300_H89up_promoter_olp), "P300_H89up_promoter_olp.csv", sep=",")
```

 - Which ATAC-seq regions overlap P300 regions?
 
```R
seqlevelsStyle(P300_H89down) <- "UCSC"
seqlevelsStyle(P300_H89up)   <- "UCSC"

P300_H89down_dyn_inc_olp <- subsetByOverlaps(P300_H89down, dyn_inc_GR)
P300_H89up_dyn_inc_olp   <- subsetByOverlaps(P300_H89up, dyn_inc_GR)
fwrite(as.data.table(P300_H89down_dyn_inc_olp), "P300_H89down_dyn_inc_olp.csv", sep=",")
fwrite(as.data.table(P300_H89up_dyn_inc_olp), "P300_H89up_dyn_inc_olp.csv", sep=",")

P300_H89down_dyn_dec_olp <- subsetByOverlaps(P300_H89down, dyn_dec_GR)
P300_H89up_dyn_dec_olp   <- subsetByOverlaps(P300_H89up, dyn_dec_GR)
fwrite(as.data.table(P300_H89down_dyn_dec_olp), "P300_H89down_dyn_dec_olp.csv", sep=",")
fwrite(as.data.table(P300_H89up_dyn_dec_olp), "P300_H89up_dyn_dec_olp.csv", sep=",")
```

 - Of the ATAC-seq regions that overlap P300 ChIP regions, which of those associated target genes are present in the marker lists from scRNA-seq?
 - Here, look at which ChIP regions overlap the RNAseq markers
 
```R
P300_H89down_trt_markers_olp <- subsetByOverlaps(
    P300_H89down, marker_files_annotated_GR[["trt_markers"]])
P300_H89up_trt_markers_olp   <- subsetByOverlaps(
    P300_H89up, marker_files_annotated_GR[["trt_markers"]])
```

 - Here, look at which ChIP-ATAC overlap regions associated target genes are present in the RNAseq marker list

```R
# Increased ATAC regions
markers_olp_P300_H89down_and_incATAC <- dplyr::intersect(marker_files_annotated_GR[["trt_markers"]]$gene, P300_H89down_dyn_inc_olp$symbol)
P300_H89down_incATAC_RNA_olp <- marker_files_annotated_GR[["trt_markers"]][marker_files_annotated_GR[["trt_markers"]]$gene %in% markers_olp_P300_H89down_and_incATAC,]
fwrite(as.data.table(P300_H89down_incATAC_RNA_olp), "P300_H89down-incATAC-RNA_olp.csv", sep=",")
# only gene here is Cacna1a. A calcium channel. P300 binding down, but increased accessibility/

markers_olp_P300_H89up_and_incATAC <- dplyr::intersect(marker_files_annotated_GR[["trt_markers"]]$gene, P300_H89up_dyn_inc_olp$symbol)
P300_H89up_incATAC_RNA_olp <- marker_files_annotated_GR[["trt_markers"]][marker_files_annotated_GR[["trt_markers"]]$gene %in% markers_olp_P300_H89up_and_incATAC,]
fwrite(as.data.table(P300_H89up_incATAC_RNA_olp), "P300_H89up-incATAC-RNA_olp.csv", sep=",")

# Decreased ATAC regions
markers_olp_P300_H89down_and_decATAC <- dplyr::intersect(marker_files_annotated_GR[["trt_markers"]]$gene, P300_H89down_dyn_dec_olp$symbol)
P300_H89down_decATAC_RNA_olp <- marker_files_annotated_GR[["trt_markers"]][marker_files_annotated_GR[["trt_markers"]]$gene %in% markers_olp_P300_H89down_and_decATAC,]
fwrite(as.data.table(P300_H89down_decATAC_RNA_olp), "P300_H89down-decATAC-RNA_olp.csv", sep=",")

markers_olp_P300_H89up_and_decATAC <- dplyr::intersect(marker_files_annotated_GR[["trt_markers"]]$gene, P300_H89up_dyn_dec_olp$symbol)
P300_H89up_decATAC_RNA_olp <- marker_files_annotated_GR[["trt_markers"]][marker_files_annotated_GR[["trt_markers"]]$gene %in% markers_olp_P300_H89up_and_decATAC,]
fwrite(as.data.table(P300_H89up_decATAC_RNA_olp), "P300_H89up-decATAC-RNA_olp.csv", sep=",")
```

Interpretation example:

"P300_H89down-incATAC-RNA_olp.csv" 
  - These are regions where P300 is reduced following H89 treatment
  - These are regions that have increased accessibility from the ATAC-seq data following H89 treatment
  - These are the genes that are marker genes that separate scRNA-seq between H89 and DMSO

 - Here, look at which ChIP-ATAC overlap regions associated target genes are present in the `Ren1_LO` cluster's RNAseq marker list

```R
# Increased ATAC regions
markers_olp_P300_H89down_and_incATAC <- dplyr::intersect(marker_files_annotated_GR[["ren1_lo_ren1_lo_trt_markers"]]$gene, P300_H89down_dyn_inc_olp$symbol)
P300_H89down_incATAC_RNA_olp <- marker_files_annotated_GR[["ren1_lo_trt_markers"]][marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene %in% markers_olp_P300_H89down_and_incATAC,]
fwrite(as.data.table(P300_H89down_incATAC_RNA_olp), "P300_H89down-incATAC-Ren1-LO-RNA_olp.csv", sep=",")

markers_olp_P300_H89up_and_incATAC <- dplyr::intersect(marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene, P300_H89up_dyn_inc_olp$symbol)
P300_H89up_incATAC_RNA_olp <- marker_files_annotated_GR[["ren1_lo_trt_markers"]][marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene %in% markers_olp_P300_H89up_and_incATAC,]
fwrite(as.data.table(P300_H89up_incATAC_RNA_olp), "P300_H89up-incATAC-Ren1-LO-RNA_olp.csv", sep=",")

# Decreased ATAC regions
markers_olp_P300_H89down_and_decATAC <- dplyr::intersect(marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene, P300_H89down_dyn_dec_olp$symbol)
P300_H89down_decATAC_RNA_olp <- marker_files_annotated_GR[["ren1_lo_trt_markers"]][marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene %in% markers_olp_P300_H89down_and_decATAC,]
fwrite(as.data.table(P300_H89down_decATAC_RNA_olp), "P300_H89down-decATAC-Ren1-LO-RNA_olp.csv", sep=",")

markers_olp_P300_H89up_and_decATAC <- dplyr::intersect(marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene, P300_H89up_dyn_dec_olp$symbol)
P300_H89up_decATAC_RNA_olp <- marker_files_annotated_GR[["ren1_lo_trt_markers"]][marker_files_annotated_GR[["ren1_lo_trt_markers"]]$gene %in% markers_olp_P300_H89up_and_decATAC,]
fwrite(as.data.table(P300_H89up_decATAC_RNA_olp), "P300_H89up-decATAC-Ren1-LO-RNA_olp.csv", sep=",")
```

Okay. So, can I identify a subset of genes that under H89 treatment are:
1. elevated expression
2. increased accessibility
3. differential P300 binding and H3K27Ac

Example:
`P300_H89up-incATAC-Ren1-LO-RNA_olp.csv` | `P300_H89up_incATAC_RNA_olp`
 - P300 binding is increased following H89 in these regions that are increased in accessibility and overlap this subset of genes that are markers in the `Ren1_LO cluster` from the scRNA-seq. 

`H3K27Ac_H89up-incATAC-Ren1-LO-RNA_olp.csv` | `H3K27Ac_H89up_incATAC_RNA_olp`
 - This represents regions with increased H3K27Ac and accessibility that overlap markers in the `Ren1_LO cluster` from the scRNA-seq.
 
 - Do these genes overlap?

```R
# All up
tmp <- merge(as.data.table(P300_H89up_incATAC_RNA_olp),
             as.data.table(H3K27Ac_H89up_incATAC_RNA_olp))
H89up_P300_H3K27Ac_incATAC_Ren1LO_genes <- intersect(P300_H89up_incATAC_RNA_olp$gene, H3K27Ac_H89up_incATAC_RNA_olp$gene)
fwrite(tmp, "H89up_P300-H3K27Ac_incATAC_Ren1-LO_olp.csv", sep=",")

# All down
tmp <- merge(as.data.table(P300_H89down_decATAC_RNA_olp),
             as.data.table(H3K27Ac_H89down_decATAC_RNA_olp))
H89up_P300_H3K27Ac_decATAC_Ren1LO_genes <- intersect(P300_H89down_decATAC_RNA_olp$gene, H3K27Ac_H89down_decATAC_RNA_olp$gene)
fwrite(tmp, "H89down_P300-H3K27Ac_decATAC_Ren1-LO_olp.csv", sep=",")
```

### 3i. Plot overlapping genes/regions/marks

```R
library(dplyr)
library(patchwork)
library(svglite)
library(ggplot2)
library(sctransform) # must load before Seurat
library(Seurat)
library(glmGamPoi)
library(data.table)
library(harmony)
library(fastSave)
library(ggstatsplot)

set.seed(99)

H89up_P300_H3K27Ac_incATAC_Ren1LO_genes <- c(
    "Btg1", "Asap1", "Il17ra", "Rsu1", "Scd1",
    "Tbc1d2", "Anxa4", "Ndrg4", "Tm4sf1")

cluster <- "Ren1_LO"
cluster_GEX_dt <- GEX_dt[named_cluster == cluster,]

for (gene in H89up_P300_H3K27Ac_incATAC_Ren1LO_genes) {
    g1 <- ggstatsplot::ggbetweenstats(
      data  = cluster_GEX_dt,
      x     = group,
      y     = !!rlang::sym(gene),
      title = gene,
      ggtheme = custom_theme() + theme(panel.grid.major = element_line(
                                        color = "gray",
                                        linewidth = 0.1,
                                        linetype = 2))
    )
    g1 + ylim(0,7)

    svglite(paste0(gene, "_SCT_ggbetweenstats_", stringr::str_replace(cluster, "\\|", "-"), ".svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(g1 + ylim(0,7))
    dev.off()
}

for (gene in unique(H89up_P300_H3K27Ac_incATAC_Ren1LO_genes)) {
    g1 <- ggstatsplot::ggbetweenstats(
      data  = GEX_dt,
      x     = group,
      y     = !!rlang::sym(gene),
      title = gene,
      ggtheme = custom_theme() + theme(panel.grid.major = element_line(
                                        color = "gray",
                                        linewidth = 0.1,
                                        linetype = 2))
    )
    g1 + ylim(0,max_GEX_val)

    svglite(paste0(gene, "_CTRL_SCT_ggbetweenstats.svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(g1 + ylim(0,max_GEX_val))
    dev.off()
}

H89down_P300_H3K27Ac_decATAC_Ren1LO_genes <- c("Pfkfb3", "Cdk1", "Pkdcc")
for (gene in unique(H89down_P300_H3K27Ac_decATAC_Ren1LO_genes)) {
    g1 <- ggstatsplot::ggbetweenstats(
      data  = GEX_dt,
      x     = group,
      y     = !!rlang::sym(gene),
      title = gene,
      ggtheme = custom_theme() + theme(panel.grid.major = element_line(
                                        color = "gray",
                                        linewidth = 0.1,
                                        linetype = 2))
    )
    g1 + ylim(0,max_GEX_val)

    svglite(paste0(gene, "_CTRL_SCT_ggbetweenstats.svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(g1 + ylim(0,max_GEX_val))
    dev.off()
}

for (gene in H89down_P300_H3K27Ac_decATAC_Ren1LO_genes) {
    g1 <- ggstatsplot::ggbetweenstats(
      data  = cluster_GEX_dt,
      x     = group,
      y     = !!rlang::sym(gene),
      title = gene,
      ggtheme = custom_theme() + theme(panel.grid.major = element_line(
                                        color = "gray",
                                        linewidth = 0.1,
                                        linetype = 2))
    )
    g1 + ylim(0,7)

    svglite(paste0(gene, "_SCT_ggbetweenstats_", stringr::str_replace(cluster, "\\|", "-"), ".svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(g1 + ylim(0,7))
    dev.off()
}
```

### 3j. Move files to public genome_browser/ locations

```console
cd $PROCESSED/H89/integration/

export BIGWIG_DIR="$PROCESSED/H89/FIMO/tomtom/fimo_composites/main_figure_beds/"
export WWW_DIR="$PROCESSED/H89/integration/genome_browser/trackHub/mm10/"

```

### 3k. Convert the BED files to bigBed
 
 - First, ensure BED files are in proper format

```R
set.seed(99)

library(data.table)
library(rtracklayer)
library(GenomeInfoDb)

integration_dir <- "$PROCESSED/H89/integration/"

chrom_sizes <- fread('mm10.chrom.sizes')
colnames(chrom_sizes) <- c("chr", "end")
chrom_sizes[, start:=1]
chrom_GR <- makeGRangesFromDataFrame(chrom_sizes, keep.extra.columns=TRUE)

bed_files <- c("dynamic_peaks_DECREASED.bed", 
               "dynamic_peaks_INCREASED.bed",
               "dynamic_peaks_ALL.bed", 
               "nondynamic_peaks.bed")

bed_file <- fread(file.path(integration_dir, bed_files[1]),
                  col.names=c("chr", "start", "end", "name", "score", "strand"))
bed_GR <- makeGRangesFromDataFrame(bed_file, keep.extra.columns=TRUE)
range(bed_GR$score)
seqlevelsStyle(bed_GR) <- "UCSC"
# now remove all ranges that are not 'within' the chromosome boundaries
bed_GR <- subsetByOverlaps(bed_GR, chrom_GR, type="within")
export.bed(sort(bed_GR), "As4.1_DMSO-v-H89_decreased-dynamic-peaks.bed")

bed_file <- fread(file.path(integration_dir, bed_files[2]),
                  col.names=c("chr", "start", "end", "name", "score", "strand"))
bed_GR <- makeGRangesFromDataFrame(bed_file, keep.extra.columns=TRUE)
range(bed_GR$score)
seqlevelsStyle(bed_GR) <- "UCSC"
# now remove all ranges that are not 'within' the chromosome boundaries
bed_GR <- subsetByOverlaps(bed_GR, chrom_GR, type="within")
export.bed(sort(bed_GR), "As4.1_DMSO-v-H89_increased-dynamic-peaks.bed")

bed_file <- fread(file.path(integration_dir, bed_files[3]),
                  col.names=c("chr", "start", "end"))
bed_GR <- makeGRangesFromDataFrame(bed_file)
seqlevelsStyle(bed_GR) <- "UCSC"
# now remove all ranges that are not 'within' the chromosome boundaries
bed_GR <- subsetByOverlaps(bed_GR, chrom_GR, type="within")
export.bed(sort(bed_GR), "As4.1_DMSO-v-H89_all-dynamic-peaks.bed")

bed_file <- fread(file.path(integration_dir, bed_files[4]),
                  col.names=c("chr", "start", "end"))
bed_GR <- makeGRangesFromDataFrame(bed_file)
seqlevelsStyle(bed_GR) <- "UCSC"
# now remove all ranges that are not 'within' the chromosome boundaries
bed_GR <- subsetByOverlaps(bed_GR, chrom_GR, type="within")
export.bed(sort(bed_GR), "As4.1_DMSO-v-H89_nondynamic-peaks.bed")

q()
```

 - Now, convert the bed files to bigBed format
 
```console
LC_COLLATE=C sort -k1,1 -k2,2n As4.1_DMSO-v-H89_all-dynamic-peaks.bed > tmp.bed && mv tmp.bed As4.1_DMSO-v-H89_all-dynamic-peaks.bed
bedToBigBed -type=bed6 As4.1_DMSO-v-H89_all-dynamic-peaks.bed mm10.chrom.sizes As4.1_DMSO-v-H89_all-dynamic-peaks.bb

LC_COLLATE=C sort -k1,1 -k2,2n As4.1_DMSO-v-H89_decreased-dynamic-peaks.bed  > tmp.bed && mv tmp.bed As4.1_DMSO-v-H89_decreased-dynamic-peaks.bed 
bedToBigBed -type=bed6 As4.1_DMSO-v-H89_decreased-dynamic-peaks.bed mm10.chrom.sizes As4.1_DMSO-v-H89_decreased-dynamic-peaks.bb

LC_COLLATE=C sort -k1,1 -k2,2n As4.1_DMSO-v-H89_increased-dynamic-peaks.bed > tmp.bed && mv tmp.bed As4.1_DMSO-v-H89_increased-dynamic-peaks.bed
bedToBigBed -type=bed6 As4.1_DMSO-v-H89_increased-dynamic-peaks.bed mm10.chrom.sizes As4.1_DMSO-v-H89_increased-dynamic-peaks.bb

LC_COLLATE=C sort -k1,1 -k2,2n As4.1_DMSO-v-H89_nondynamic-peaks.bed > tmp.bed && mv tmp.bed As4.1_DMSO-v-H89_nondynamic-peaks.bed
bedToBigBed -type=bed6 As4.1_DMSO-v-H89_nondynamic-peaks.bed mm10.chrom.sizes As4.1_DMSO-v-H89_nondynamic-peaks.bb
```

Copy bigBed files to /www/
```console
export WWW_DIR="$PROCESSED/H89/integration/genome_browser/trackHub/mm10/"

cp As4.1_DMSO-v-H89_decreased-dynamic-peaks.bb ${WWW_DIR}
cp As4.1_DMSO-v-H89_increased-dynamic-peaks.bb ${WWW_DIR}
cp As4.1_DMSO-v-H89_all-dynamic-peaks.bb ${WWW_DIR}
cp As4.1_DMSO-v-H89_nondynamic-peaks.bb ${WWW_DIR}
```

Copy the ATAC signal tracks to the genome_browser/ location
```console
cp $PROCESSED/H89/bam_files/As4.1_H89_merged.bw ${WWW_DIR}
```
