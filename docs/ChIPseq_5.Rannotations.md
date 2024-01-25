### 5. Check whether ChIP-seq regions are preferentially enriched in promoters

#### 5a. H3K27Ac

```R
project_name <- "H3K27Ac"

txdb      <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)

seqlevelsStyle(promoters) <- "Ensembl"

dba_H3K27Ac_peaks_DMSOenriched <- dba_H3K27Ac_peaksAnno[dba_H3K27Ac_peaksAnno$Fold > 0,]
dba_H3K27Ac_peaks_H89enriched  <- dba_H3K27Ac_peaksAnno[dba_H3K27Ac_peaksAnno$Fold < 0,]

export.bed(dba_H3K27Ac_peaks_H89enriched,
           con=paste0(project_name, "_DB/", project_name,
               "_H89_enriched_peaks.bed"))

export.bed(dba_H3K27Ac_peaks_DMSOenriched,
           con=paste0(project_name, "_DB/", project_name,
               "_DMSO_enriched_peaks.bed"))

prom_overlap <- findOverlaps(dba_H3K27Ac_peaks_H89enriched, promoters)

cat(sprintf("%d of %d promoters are overlapped by an enriched region.",
   length(unique(subjectHits(prom_overlap))), length(promoters)))
```
3567 of 24528 promoters are overlapped by an enriched region.

 - Now the inverse, which H3K27ac enriched regions overlap promoters?
```R
prom_overlap2 <- findOverlaps(promoters, dba_H3K27Ac_peaks_H89enriched)

cat(sprintf( "%d of %d enriched regions overlap a promoter.",
   length( unique( subjectHits(prom_overlap2))), length(dba_H3K27Ac_peaks_H89enriched)))
```
3025 of 11040 enriched regions overlap a promoter.

```R
promotersDT <- as.data.table(promoters)
promoter_length_by_chr <- promotersDT[ ,list(sum=sum(width)), by=seqnames]

tmp <- copy(genome)
seqlevelsStyle(tmp) <- "Ensembl"
chr_lengths <- as.data.table(seqlengths(tmp), keep.rownames=TRUE)
colnames(chr_lengths) <- c("seqnames", "sum")
rm(tmp)

promoter_seqlength <- merge(chr_lengths, promoter_length_by_chr, by="seqnames")

# Which fraction of the chromsome is this?
promoter_seqlength$promoter_fraction <- promoter_seqlength$sum.y / promoter_seqlength$sum.x

# Keep only chromosomes present in the dba_object
promoter_seqlength <- promoter_seqlength[promoter_seqlength$seqnames %in% 
    seqnames(dba_H3K27Ac_peaks_H89enriched),]

# Is this significantly enriched? Yes
binom.test(length(unique(subjectHits(prom_overlap2))), length(dba_H3K27Ac_peaks_H89enriched), mean(promoter_seqlength$promoter_fraction))
```
        Exact binomial test

data:  length(unique(subjectHits(prom_overlap2))) and length(dba_H3K27Ac_peaks_H89enriched)
number of successes = 3025, number of trials = 11040, p-value < 2.2e-16
alternative hypothesis: true probability of success is not equal to 0.01886861
95 percent confidence interval:
 0.2656991 0.2824275
sample estimates:
probability of success
             0.2740036

#### 5b. P300

```R
project_name <- "P300"

txdb      <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)

seqlevelsStyle(promoters) <- "Ensembl"

dba_P300_peaks_DMSOenriched <- dba_P300_peaksAnno[dba_P300_peaksAnno$Fold > 0,]
dba_P300_peaks_H89enriched  <- dba_P300_peaksAnno[dba_P300_peaksAnno$Fold < 0,]

export.bed(dba_P300_peaks_H89enriched,
           con=paste0(project_name, "_DB/", project_name,
               "_H89_enriched_peaks.bed"))

export.bed(dba_P300_peaks_DMSOenriched,
           con=paste0(project_name, "_DB/", project_name,
               "_DMSO_enriched_peaks.bed"))


prom_overlap <- findOverlaps(dba_P300_peaks_H89enriched, promoters)

cat(sprintf("%d of %d promoters are overlapped by an enriched region.",
   length(unique(subjectHits(prom_overlap))), length(promoters)))
```
17 of 24528 promoters are overlapped by an enriched region.

 - Now the inverse, which P300 enriched regions overlap promoters?
```R
prom_overlap2 <- findOverlaps(promoters, dba_P300_peaks_H89enriched)

cat(sprintf( "%d of %d enriched regions overlap a promoter.",
   length( unique( subjectHits(prom_overlap2))), length(dba_P300_peaks_H89enriched)))
```
17 of 1227 enriched regions overlap a promoter.

```R
promotersDT <- as.data.table(promoters)
promoter_length_by_chr <- promotersDT[ ,list(sum=sum(width)), by=seqnames]

tmp <- copy(genome)
seqlevelsStyle(tmp) <- "Ensembl"
chr_lengths <- as.data.table(seqlengths(tmp), keep.rownames=TRUE)
colnames(chr_lengths) <- c("seqnames", "sum")
rm(tmp)

promoter_seqlength <- merge(chr_lengths, promoter_length_by_chr, by="seqnames")

# Which fraction of the chromsome is this?
promoter_seqlength$promoter_fraction <- promoter_seqlength$sum.y / promoter_seqlength$sum.x

# Keep only chromosomes present in the dba_object
promoter_seqlength <- promoter_seqlength[promoter_seqlength$seqnames %in% 
    seqnames(dba_P300_peaks_H89enriched),]

# Is this significantly enriched? NO
binom.test(length(unique(subjectHits(prom_overlap2))), length(dba_P300_peaks_H89enriched), mean(promoter_seqlength$promoter_fraction))
```
        Exact binomial test

data:  length(unique(subjectHits(prom_overlap2))) and length(dba_P300_peaks_H89enriched)
number of successes = 17, number of trials = 1227, p-value = 0.2466
alternative hypothesis: true probability of success is not equal to 0.01886861
95 percent confidence interval:
 0.008091166 0.022090670
sample estimates:
probability of success
            0.01385493