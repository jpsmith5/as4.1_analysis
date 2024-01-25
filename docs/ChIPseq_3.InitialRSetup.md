## 3. Analyze outputs in R

### 3a. Create samplesheet

```console
tr '\t' ',' < $PROCESSED/nf-core/chipseq/diffbind_H89.csv > tmp.csv && mv tmp.csv $PROCESSED/nf-core/chipseq/diffbind_H89.csv
```

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
library(scales)
library(stringi)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)
library(rtracklayer)

txdb   <- TxDb.Mmusculus.UCSC.mm10.knownGene
genome <- BSgenome.Mmusculus.UCSC.mm10

samples <- read.csv("$PROCESSED/nf-core/chipseq/diffbind_H89.csv")
```

### 3b. Perform analysis

```R
samples_H3K27Ac <- samples[1:4,]
samples_P300 <- samples[5:8,]

default_bam_dir <- "$PROCESSED/nf-core/chipseq/H89/bwa/mergedLibrary/"
chr_bam_dir     <- "$PROCESSED/nf-core/chipseq/bamfiles/"
chr_ctl_bam_dir <- "$PROCESSED/nf-core/chipseq/bamfiles/control/"
samples_H3K27Ac$bamReads <- str_replace(samples_H3K27Ac$bamReads,
                                        default_bam_dir,
                                        chr_bam_dir)
samples_H3K27Ac$bamControl <- str_replace(samples_H3K27Ac$bamControl,
                                          default_bam_dir,
                                          chr_ctl_bam_dir)
samples_H3K27Ac$bamReads <- str_replace(samples_H3K27Ac$bamReads,
                                        "\\.bam",
                                        "_chr.bam")
samples_H3K27Ac$bamControl <- str_replace(samples_H3K27Ac$bamControl,
                                          "\\.bam",
                                          "_chr.bam")
```

### 3c. Create some easy to use group functions for DiffBind analysis.

```R
calcDBAobject <- function(samples, project_name) {
    dir.create(file.path(getwd(), paste0(project_name, "_DB/")),
                         showWarnings = FALSE)

    dba_object <- dba(sampleSheet=samples)
    # Convert to summarizedExperiment object to see the chromosome naming style
    sset <- dba(dba_object,bSummarizedExperiment=TRUE)
    seq_style <- seqlevelsStyle(sset)[[length(seqlevelsStyle(sset))]]

    # Load the blacklist object `mm10.blacklist`
    load(system.file("data/mm10.blacklist.RData",package="GreyListChIP"),
         envir = environment())
    # Convert to Ensembl to match input when using the base input files,
    # not the swapped alternatives
    seqlevelsStyle(mm10.blacklist) <- seq_style

    dba_object <- dba.blacklist(dba_object,
                                blacklist=mm10.blacklist,
                                greylist=FALSE)
    dba_object <- dba.count(dba_object)
    dba_object <- dba.normalize(dba_object, method=DBA_ALL_METHODS)
    dba_object <- dba.contrast(dba_object,
                               categories=DBA_CONDITION,
                               minMembers=2,
                               reorderMeta = list(Condition=c("H89")))
    dba_object <- dba.analyze(dba_object,
                              method=DBA_ALL_METHODS,
                              bParallel = TRUE)
    info <- dba.show(dba_object)

    libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,
                      PeakReads=round(info$Reads * info$FRiP))
    rownames(libsizes) <- info$ID

    norm_object <- dba.normalize(dba_object, bRetrieve=TRUE)
    normlibs <- cbind(FullLibSize=norm_object$lib.sizes,
                      NormFacs=norm_object$norm.factors,
                      NormLibSize=round(
                        norm_object$lib.sizes/norm_object$norm.factors))
    rownames(normlibs) <- info$ID

    rm(consensus_peaks)
    for (i in 1:4) {
        name  <- str_replace(dba_object$samples$SampleID[[i]], "As4.1_", "")
        peaks <- dba_object$peaks[[i]]
        colnames(peaks) <- c("Chr", "Start", "End", paste0(name, "_Score"),
                             paste0(name, "_RPKM"), paste0(name, "_Reads"),
                             paste0(name, "_cRPKM"), paste0(name, "_cReads"))
        GR <- makeGRangesFromDataFrame(peaks, keep.extra.columns=TRUE)
        if (exists("consensus_peaks")) {
            consensus_peaks <- merge(consensus_peaks, GR)
        } else {
            consensus_peaks <- GR
        }
    }

    export.bed(consensus_peaks,
        con=paste0(project_name, "_DB/", project_name, "_consensus_peaks.bed")
    )

    return(dba_object)
}

plotDBAobject <- function(dba_object, project_name) {
    svglite(paste0(project_name, "_DB/", project_name,
                   "_H89-v-DMSO_DBplot.svg"))
    plot(dba_object, contrast=1)
    dev.off()

    svglite(paste0(project_name, "_DB/", project_name,
                   "_H89-v-DMSO_VennPlot.svg"))
    dba.plotVenn(dba_object, contrast=1, bDB=TRUE,
                       bGain=TRUE, bLoss=TRUE, bAll=FALSE)
    dev.off()

    svglite(paste0(project_name, "_DB/", project_name,
                   "_H89-v-DMSO_PCAplot.svg"))
    dba.plotPCA(dba_object, attributes=c(DBA_FACTOR,DBA_CONDITION))
    dev.off()

    svglite(paste0(project_name, "_DB/", project_name,
                   "_H89-v-DMSO_MAplot.svg"))
    dba.plotMA(dba_object)
    dev.off()

    svglite(paste0(project_name, "_DB/", project_name,
                   "_H89-v-DMSO_VolcanoPlot.svg"))
    dba.plotVolcano(dba_object)
    dev.off()

    # boxplots to view how read distributions differ between classes of binding sites
    svglite(paste0(project_name, "_DB/", project_name,
                   "_H89-v-DMSO_BoxPlot.svg"))
    pvals <- dba.plotBox(dba_object, contrast=1) # Only 1 contrast anyways...
    dev.off()

    corvals <- dba.plotHeatmap(dba_object)
    hmap    <- colorRampPalette(c("blue", "white", "yellow"))(n = 13)
    svglite(paste0(project_name, "_DB/", project_name,
                   "_H89-v-DMSO_Heatmap.svg"))
    readscores <- dba.plotHeatmap(dba_object, contrast=1, correlations=FALSE,
                                  scale="row", colScheme = hmap)
    dev.off()

    profiles <- dba.plotProfile(dba_object)
    svglite(paste0(project_name, "_DB/", project_name,
                   "_H89-v-DMSO_ProfilesPlot.svg"))
    dba.plotProfile(profiles)
    dev.off()

    olap.rate <- dba.overlap(dba_object,mode=DBA_OLAP_RATE)
    plot(olap.rate,type='b',ylab='# peaks',
         xlab='Overlap at least this many peaksets')
    names(dba_object$masks)
    dba.overlap(dba_object,dba_object$masks$H89, mode=DBA_OLAP_RATE)
    dba.plotVenn(dba_object, dba_object$masks$H89)

    svglite(paste0(project_name, "_DB/", project_name,
                   "_H89-v-DMSO_ConsensusPeakVenn.svg"))
    dba.plotVenn(dba_object, dba_object$masks$Consensus)
    dev.off()
}
```