### 9. Now, only keep ChIP-seq peaks that overlap the ATAC-seq peaks

#### 9a. H3K27Ac

```R
project_name <- "H3K27Ac"

# Which regions in the ChIP data overlap ATAC-seq regions
true_peak_olp <- findOverlaps(H89_ATACseq_consensusGR,
                              dba_H3K27Ac_peaks_H89enriched)
# Subset ChIP data on these hits
dba_H3K27Ac_H89_true_peaks <- dba_H3K27Ac_peaks_H89enriched[
    subjectHits(true_peak_olp),]

true_peaksDT <- as.data.table(dba_H3K27Ac_H89_true_peaks)
true_peaksDT <- dplyr::distinct(true_peaksDT, .keep_all = TRUE)
true_peaks_H3K27Ac_H89 <- makeGRangesFromDataFrame(true_peaksDT,
    keep.extra.columns=TRUE)

true_peaks_H3K27Ac_H89[grep("Ren1", true_peaks_H3K27Ac_H89$symbol),]
# N/A
true_peaks_H3K27Ac_H89[grep("Nfix", true_peaks_H3K27Ac_H89$symbol),]
# present

export.bed(true_peaks_H3K27Ac_H89,
           con=paste0(project_name, "_DB/", project_name,
               "_H89_true_peaks.bed"))
fwrite(as.data.table(true_peaks_H3K27Ac_H89),
       paste0(project_name, "_DB/", project_name,
              "_H89_true_peaks.csv"))

# Which regions in the DMSO ChIP data overlap DMSO ATAC-seq regions
true_peak_olp <- findOverlaps(DMSO_ATACseq_consensusGR,
                              dba_H3K27Ac_peaks_DMSOenriched)
# Subset ChIP data on these hits
dba_H3K27Ac_DMSO_true_peaks <- dba_H3K27Ac_peaks_DMSOenriched[subjectHits(true_peak_olp),]

true_peaksDT <- as.data.table(dba_H3K27Ac_DMSO_true_peaks)
true_peaksDT <- dplyr::distinct(true_peaksDT, .keep_all = TRUE)
true_peaks_H3K27Ac_DMSO <- makeGRangesFromDataFrame(true_peaksDT, keep.extra.columns=TRUE)

true_peaks_H3K27Ac_DMSO[grep("Ren1", true_peaks_H3K27Ac_DMSO$symbol),]
# present
true_peaks_H3K27Ac_DMSO[grep("Nfix", true_peaks_H3K27Ac_DMSO$symbol),]
# N/A

export.bed(true_peaks_H3K27Ac_DMSO,
           con=paste0(project_name, "_DB/", project_name,
               "_DMSO_true_peaks.bed"))
fwrite(as.data.table(true_peaks_H3K27Ac_DMSO),
       paste0(project_name, "_DB/", project_name,
              "_DMSO_true_peaks.csv"))
```

#### 9b. P300

```R
project_name <- "P300"

# Which regions in the ChIP data overlap ATAC-seq regions
true_peak_olp <- findOverlaps(H89_ATACseq_consensusGR,
                              dba_P300_peaks_H89enriched)
# Subset ChIP data on these hits
dba_P300_H89_true_peaks <- dba_P300_peaks_H89enriched[
    subjectHits(true_peak_olp),]

true_peaksDT <- as.data.table(dba_P300_H89_true_peaks)
true_peaksDT <- dplyr::distinct(true_peaksDT, .keep_all = TRUE)
true_peaks_P300_H89 <- makeGRangesFromDataFrame(true_peaksDT,
    keep.extra.columns=TRUE)

true_peaks_P300_H89[grep("Ren1", true_peaks_P300_H89$symbol),]
# N/A
true_peaks_P300_H89[grep("Nfix", true_peaks_P300_H89$symbol),]
# present

export.bed(true_peaks_P300_H89,
           con=paste0(project_name, "_DB/", project_name,
               "_H89_true_peaks.bed"))
fwrite(as.data.table(true_peaks_P300_H89),
       paste0(project_name, "_DB/", project_name,
              "_H89_true_peaks.csv"))

# Which regions in the DMSO ChIP data overlap DMSO ATAC-seq regions
true_peak_olp <- findOverlaps(DMSO_ATACseq_consensusGR,
                              dba_P300_peaks_DMSOenriched)
# Subset ChIP data on these hits
dba_P300_DMSO_true_peaks <- dba_P300_peaks_DMSOenriched[subjectHits(true_peak_olp),]

true_peaksDT <- as.data.table(dba_P300_DMSO_true_peaks)
true_peaksDT <- dplyr::distinct(true_peaksDT, .keep_all = TRUE)
true_peaks_P300_DMSO <- makeGRangesFromDataFrame(true_peaksDT, keep.extra.columns=TRUE)

true_peaks_P300_DMSO[grep("Ren1", true_peaks_P300_DMSO$symbol),]
# present
true_peaks_P300_DMSO[grep("Nfix", true_peaks_P300_DMSO$symbol),]
# N/A

export.bed(true_peaks_P300_DMSO,
           con=paste0(project_name, "_DB/", project_name,
               "_DMSO_true_peaks.bed"))
fwrite(as.data.table(true_peaks_P300_DMSO),
       paste0(project_name, "_DB/", project_name,
              "_DMSO_true_peaks.csv"))
```