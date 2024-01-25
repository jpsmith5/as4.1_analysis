### 8a. Subset ChIP-seq peaks by overlapping to ATAC-seq open chromatin regions

Look at regions of the ChIPseq data that overlap with ATACseq.

```console
cd $PROCESSED/H89/
mkdir peak_files/
find . -type f -name '*normalized.narrowPeak' -exec cp {} peak_files/ \;
```

### 8b. Load ATAC-seq peak sets

```R
H89_ATACseq_path  <- "$PROCESSED/H89/peak_files/"
H89_ATACseq_files <- dir(H89_ATACseq_path, "narrowPeak")
H89_ATACseq_data  <- lapply(file.path(H89_ATACseq_path, H89_ATACseq_files),
                            fread)
names(H89_ATACseq_data) <- gsub(".narrowPeak", "", H89_ATACseq_files)
for (sample_name in names(H89_ATACseq_data)) {
    colnames(H89_ATACseq_data[[sample_name]]) <- c(
        "chr", "chromStart", "chromEnd", "name", "score",
        "strand", "signalValue", "pValue", "qValue", "peak")
    H89_ATACseq_data[[sample_name]]$score <- 
        rescale(log(H89_ATACseq_data[[sample_name]]$score), to = c(0, 1000))
    H89_ATACseq_data[[sample_name]] <- 
        makeGRangesFromDataFrame(H89_ATACseq_data[[sample_name]],
                                 keep.extra.columns=TRUE)
}
```

### 8c. Combine treatment groups in consensus peak sets

```R
collapsePeaks <- function(peak_files, sample_names, chr_sizes,
						  min_samples=2, min_score=5, min_olap=1) {
    # create combined peaks
    peaks <- rbindlist(lapply(peak_files, fread), idcol="file")
    if (ncol(peaks) == 7) {
        colnames(peaks) <- c("file", "chr", "start", "end",
                             "name", "score", "strand")
    } else if (ncol(peaks) == 11) {
        colnames(peaks) <- c("file", "chr", "start", "end",
                             "name", "score", "strand",
                             "signalValue", "pValue", "qValue", "peak")
    } else {
        warning(paste0("Peak files did not contain a recognizable number", 
                       " of columns (", ncol(peaks), ")"))
        rm(peaks)
        final <- data.table(chr=character(),
                            start=integer(),
                            end=integer(),
                            name=character(),
                            score=numeric(),
                            strand=character(),
                            signalValue=numeric(),
                            pValue=numeric(),
                            qValue=numeric(),
                            peak=integer())
        return(final)
    }
    setkey(peaks, chr, start, end)
    # keep highest scored peaks
    # split by chromosome to minimize memory requirements
    peaks_by_chr   <- split(peaks, peaks$chr)
    hit_aggregator <- function(x) {
        #message(paste0("x: ", unique(x$chr)))  # DEBUG
        peaksGR <- makeGRangesFromDataFrame(x, keep.extra.columns=FALSE)
        hitsGR  <- suppressWarnings(
            findOverlaps(peaksGR, peaksGR,
						 ignore.strand=TRUE, minoverlap=min_olap))
        hits    <- data.table::data.table(xid=queryHits(hitsGR),
                                          yid=subjectHits(hitsGR))
        setkey(hits, xid)
        scores  <- data.table(index=rep(1:nrow(x)), score=x$score)
        setkey(scores, index)
        out     <- hits[scores, nomatch=0]
        keep    <- out[out[,.I[which.max(score)],by=yid]$V1]
        indices <- unique(keep$xid)
        reduced <- x[indices,]
        reduced[start < 0, start := 0]
        return(reduced)
    }
    final <- rbindlist(lapply(peaks_by_chr, hit_aggregator))

    # can't extend past chromosome
    for (i in nrow(chr_sizes)) {
        final[chr == chr_sizes$chr[i] & end > chr_sizes$size[i],
              end := chr_sizes$size[i]]
    }

    # identify reproducible peaks
    peaks[,group := sample_names[file]]
    peaks[,file:=NULL]
    final[,file:=NULL]
    peak_list <- splitDataTable(peaks, "group")
    rm(peaks)
    invisible(gc())
    invisible(sapply(peak_list, countReproduciblePeaks, peak_DT=final))

    # keep peaks present in 2 or more individual peak sets
    # keep peaks with score per million >= 5
    final <- final[count >= min_samples & score >= min_score,]
    final[,count := NULL]
    return(final)
}

splitDataTable <- function(DT, split_factor) {
	if (is.numeric(split_factor)) {
		split_factor = colnames(DT)[split_factor]
		message("Integer split_factor, changed to: ", split_factor)
	}
	lapply( split(1:nrow(DT), DT[, get(split_factor)]), function(x) DT[x])
}

countReproduciblePeaks <- function(peak_list, peak_DT) {
    setkey(peak_DT, chr, start, end)
    setkey(peak_list, chr, start, end)
    hits <- foverlaps(peak_list, peak_DT,
                      by.x=c("chr", "start", "end"),
                      type="any", which=TRUE, nomatch=0)
    # track the number of overlaps of final peak set peaks
    if (!"count" %in% colnames(peak_DT)) {
        peak_DT[hits$yid, count := 1]
        peak_DT[is.na(get("count")), ("count") := 0]
    } else {
        peak_DT[hits$yid, count := get("count") + 1] 
    }
}

colnames(chr_lengths) <- c("chr", "size")

DMSO_ATACseq_consensus <- collapsePeaks(
    peak_files=file.path(H89_ATACseq_path, H89_ATACseq_files)[1:3],
    chr_sizes=chr_lengths,
    sample_names=str_replace(names(H89_ATACseq_data), "_peaks_normalized", ""))
DMSO_ATACseq_consensusGR <- makeGRangesFromDataFrame(DMSO_ATACseq_consensus,
    keep.extra.columns=TRUE)

H89_ATACseq_consensus <- collapsePeaks(
    peak_files=file.path(H89_ATACseq_path, H89_ATACseq_files)[4:5],
    chr_sizes=chr_lengths,
    sample_names=str_replace(names(H89_ATACseq_data), "_peaks_normalized", ""))
H89_ATACseq_consensusGR <- makeGRangesFromDataFrame(H89_ATACseq_consensus,
    keep.extra.columns=TRUE)

seqlevelsStyle(H89_ATACseq_consensusGR)  <- "Ensembl"
seqlevelsStyle(DMSO_ATACseq_consensusGR) <- "Ensembl"
```