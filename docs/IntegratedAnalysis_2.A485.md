# Integrated Analyses

## 2. A485 treated cells

Combine the A485 results with the scRNA-seq results to identify specific subsets of genomic regions and genes most affected by A485 treatment. Use the scRNA-seq CTRL (untreated) samples and evaluate the expression of the genes targeted by A485.

### 2a. Copy files for A485 comparisons

1. Copy over RNAseq files
```console
cd $PROCESSED/A485/integration_analysis/

cp $PROCESSED/scRNA/R_analysis/As4.1_H89_TRT_markers.csv .
cp $PROCESSED/scRNA/R_analysis/As4.1_H89_CTRL_markers.csv .
cp $PROCESSED/scRNA/R_analysis/markers/*.csv .
```

2. Copy over ATACseq files
```console
export ATACDIR="$PROCESSED/A485/bam_files/FIMO/tomtom/fimo_composites/"
export BAMDIR="$PROCESSED/A485/bam_files/"

cp ${ATACDIR}/cluster_bed_cluster1.bed6 ./dynamic_peaks_INCREASED.bed
cp ${ATACDIR}/cluster_bed_cluster2.bed6 ./dynamic_peaks_DECREASED.bed
cp ${BAMDIR}/dynamic_peaks.bed ./dynamic_peaks_ALL.bed
cp ${BAMDIR}/nondynamic_peaks.bed ./nondymamic_peaks.bed
cp ${BAMDIR}/all_peaks.bed ./all_peaks.bed
```

### 2b. Load data

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

directory  <- "$PROCESSED/A485/bam_files/"
tomtom_dir <- "$PROCESSED/A485/bam_files/FIMO/tomtom/"
integration_dir <- "$PROCESSED/A485/integration_analysis/"

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

### 2c. Annotate 

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

### 2d. Link ATAC-seq to scRNA-seq data 

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
fwrite(inc_REs, "A485_increased-accessibilty-regions_linked_scRNAseq.csv", sep=",")

gene_RE_table <- as.data.table(sort(table(inc_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "A485_increased-accessibilty-regions_linked_scRNAseq_table.csv", sep=",")
```

### 2e. Perform GSEA on these regions/genes

```R
GO_dir <- "$PROCESSED/A485/integration_analysis/GSEA_MSigDB/"
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
    svglite(paste0("A485_increased-accessibilty-regions_linked_scRNAseq_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}
```

### 2f. RNA marker list versus decreased dynamic ATAC-seq peaks

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
fwrite(dec_REs, "A485_decreased-accessibilty-regions_linked_scRNAseq.csv", sep=",")
```

Which genes have the most REs?

```R
gene_RE_table <- as.data.table(sort(table(dec_REs$gene_symbol)))
colnames(gene_RE_table) <- c("gene_symbol", "N")
fwrite(gene_RE_table, "A485_decreased-accessibilty-regions_linked_scRNAseq_table.csv", sep=",")
gene_RE_table[grep('Nfi', gene_RE_table$gene_symbol),]
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
    svglite(paste0("A485_decreased-accessibilty-regions_linked_scRNAseq_", GO_file, ".svg"), 
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
fwrite(ovp_genes, "A485_decreased-accessibilty-regions_linked_scRNAseq_IN_increased-accessibility-gene-list.csv", sep=",")
```

| TF_Family | bigWig_file                 |
|-----------|-----------------------------|
| AP-1      | Fra2_mm10_instances.bigWig  |
| KLF       | Klf4_mm10_instances.bigWig  |
| MAZ       | Maz_mm10_instances.bigWig   |
| SP        | Sp2_mm10_instances.bigWig   |
| P53       | Trp63_mm10_instances.bigWig |
| RUNX      | RUNX_mm10_instances.bigWig  |
| TEAD      | TEAD4_mm10_instances.bigWig |

### 2g. Move files to public genome_browser/ locations

```console
cd $PROCESSED/A485/integration/

export BIGWIG_DIR="$PROCESSED/A485/FIMO/tomtom/fimo_composites/main_figure_beds/"
export WWW_DIR="$PROCESSED/A485/integration/genome_browser/trackHub/mm10/"

cp ${BIGWIG_DIR}/Fra2_mm10_instances.bigWig ${WWW_DIR}/A485_AP1_mm10.bigWig
cp ${BIGWIG_DIR}/Klf4_mm10_instances.bigWig ${WWW_DIR}/A485_KLF_mm10.bigWig
cp ${BIGWIG_DIR}/Maz_mm10_instances.bigWig ${WWW_DIR}/A485_MAZ_mm10.bigWig
cp ${BIGWIG_DIR}/Sp2_mm10_instances.bigWig ${WWW_DIR}/A485_SP_mm10.bigWig
cp ${BIGWIG_DIR}/Trp63_mm10_instances.bigWig ${WWW_DIR}/A485_P53_mm10.bigWig
cp ${BIGWIG_DIR}/RUNX_mm10_instances.bigWig ${WWW_DIR}/A485_RUNX_mm10.bigWig
cp ${BIGWIG_DIR}/TEAD4_mm10_instances.bigWig ${WWW_DIR}/A485_TEAD_mm10.bigWig
```

### 2h. Convert the BED files to bigBed
 
 - First, ensure BED files are in proper format

```R
set.seed(99)

library(data.table)
library(rtracklayer)
library(GenomeInfoDb)

integration_dir <- "$PROCESSED/A485/integration/"

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
export.bed(sort(bed_GR), "As4.1_DMSO-v-A485_decreased-dynamic-peaks.bed")

bed_file <- fread(file.path(integration_dir, bed_files[2]),
                  col.names=c("chr", "start", "end", "name", "score", "strand"))
bed_GR <- makeGRangesFromDataFrame(bed_file, keep.extra.columns=TRUE)
range(bed_GR$score)
seqlevelsStyle(bed_GR) <- "UCSC"
# now remove all ranges that are not 'within' the chromosome boundaries
bed_GR <- subsetByOverlaps(bed_GR, chrom_GR, type="within")
export.bed(sort(bed_GR), "As4.1_DMSO-v-A485_increased-dynamic-peaks.bed")

bed_file <- fread(file.path(integration_dir, bed_files[3]),
                  col.names=c("chr", "start", "end"))
bed_GR <- makeGRangesFromDataFrame(bed_file)
seqlevelsStyle(bed_GR) <- "UCSC"
# now remove all ranges that are not 'within' the chromosome boundaries
bed_GR <- subsetByOverlaps(bed_GR, chrom_GR, type="within")
export.bed(sort(bed_GR), "As4.1_DMSO-v-A485_all-dynamic-peaks.bed")

bed_file <- fread(file.path(integration_dir, bed_files[4]),
                  col.names=c("chr", "start", "end"))
bed_GR <- makeGRangesFromDataFrame(bed_file)
seqlevelsStyle(bed_GR) <- "UCSC"
# now remove all ranges that are not 'within' the chromosome boundaries
bed_GR <- subsetByOverlaps(bed_GR, chrom_GR, type="within")
export.bed(sort(bed_GR), "As4.1_DMSO-v-A485_nondynamic-peaks.bed")

q()
```

 - Now, convert the bed files to bigBed format
 
```console
LC_COLLATE=C sort -k1,1 -k2,2n As4.1_DMSO-v-A485_all-dynamic-peaks.bed > tmp.bed && mv tmp.bed As4.1_DMSO-v-A485_all-dynamic-peaks.bed
bedToBigBed -type=bed6 As4.1_DMSO-v-A485_all-dynamic-peaks.bed mm10.chrom.sizes As4.1_DMSO-v-A485_all-dynamic-peaks.bb

LC_COLLATE=C sort -k1,1 -k2,2n As4.1_DMSO-v-A485_decreased-dynamic-peaks.bed  > tmp.bed && mv tmp.bed As4.1_DMSO-v-A485_decreased-dynamic-peaks.bed 
bedToBigBed -type=bed6 As4.1_DMSO-v-A485_decreased-dynamic-peaks.bed mm10.chrom.sizes As4.1_DMSO-v-A485_decreased-dynamic-peaks.bb

LC_COLLATE=C sort -k1,1 -k2,2n As4.1_DMSO-v-A485_increased-dynamic-peaks.bed > tmp.bed && mv tmp.bed As4.1_DMSO-v-A485_increased-dynamic-peaks.bed
bedToBigBed -type=bed6 As4.1_DMSO-v-A485_increased-dynamic-peaks.bed mm10.chrom.sizes As4.1_DMSO-v-A485_increased-dynamic-peaks.bb

LC_COLLATE=C sort -k1,1 -k2,2n As4.1_DMSO-v-A485_nondynamic-peaks.bed > tmp.bed && mv tmp.bed As4.1_DMSO-v-A485_nondynamic-peaks.bed
bedToBigBed -type=bed6 As4.1_DMSO-v-A485_nondynamic-peaks.bed mm10.chrom.sizes As4.1_DMSO-v-A485_nondynamic-peaks.bb
```

### 2i. Copy files to /www/
```console
export WWW_DIR="$PROCESSED/A485/integration/genome_browser/trackHub/mm10/"

cp As4.1_DMSO-v-A485_decreased-dynamic-peaks.bb ${WWW_DIR}
cp As4.1_DMSO-v-A485_increased-dynamic-peaks.bb ${WWW_DIR}
cp As4.1_DMSO-v-A485_all-dynamic-peaks.bb ${WWW_DIR}
cp As4.1_DMSO-v-A485_nondynamic-peaks.bb ${WWW_DIR}
```

Copy the ATAC signal tracks to the genome_browser/ location
```console
cp $PROCESSED/A485/bam_files/As4.1_A485_merged.bw ${WWW_DIR}
```

