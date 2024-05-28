# View the Ren1 locus across mouse tissues

## Obtaining public data

From [An ATAC-seq atlas of chromatin accessibility in mouse tissues](https://doi.org/10.1038/s41597-019-0071-0)
 - [Data](https://figshare.com/collections/An_ATAC-seq_atlas_of_chromatin_accessibility_in_mouse_tissues/4436264/1)

Set up R environment.

```R
library(data.table)
library(stringr)
library(rtracklayer)
```

## Load data and generate bigWig files.

```R
count_matrix <- fread('chromatin.accessibility.raw.count.txt')
count_matrix$Peak_ID <- str_replace_all(count_matrix$Peak_ID, "_random", "-random")
count_matrix[, c("chr", "start", "end") := tstrsplit(Peak_ID , "_", fixed=TRUE)]

tmat <- t(count_matrix[,2:(ncol(count_matrix)-3)])
colnames(tmat) <- count_matrix$Peak_ID
tmat <- as.data.table(tmat)
tmat[, sample:=colnames(count_matrix[,2:(ncol(count_matrix)-3)])]

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
seqlevels(txdb)

trim_dt = function(x, output){
  if(as.numeric(bg$end[i]) > seqlengths(txdb)[bg$chr[i]]) {
    bg$end[i] <- seqlengths(txdb)[bg$chr[i]]
  }
}

for (sample_name in tmat$sample) {
  bed <- as.data.table(t(tmat[sample == sample_name,]), keep.rownames=T)
  bed[, c("chr", "start", "end") := tstrsplit(rn , "_", fixed=TRUE)]
  bed <- bed[!grep("random|sample", bed$rn),]
  bed[,rn:=NULL]
  bg <- data.table(chr=bed$chr,
                   start=bed$start,
                   end=bed$end,
                   name=paste0(sample_name, "_", rep(1:nrow(bed))),
                   score=as.numeric(bed$V1),
                   strand="*")
  for(i in 1:nrow(bg)) {
    if(as.numeric(bg$end[i]) > seqlengths(txdb)[bg$chr[i]]) {
      bg$end[i] <- seqlengths(txdb)[bg$chr[i]]
    }
    if(as.numeric(bg$start[i]) > seqlengths(txdb)[bg$chr[i]]) {
      bg$start[i] <- seqlengths(txdb)[bg$chr[i]]-1
    }
  }
  bgGR <- makeGRangesFromDataFrame(bg, keep.extra.columns=T)
  genome(bgGR) <- "mm10"
  seqlevelsStyle(bgGR) <- "UCSC"
  seqlevels(bgGR) <- seqlevels(txdb)
  seqlengths(bgGR) <- seqlengths(txdb)
  bgGR <- trim(bgGR)
  export.bedGraph(bg, paste0(sample_name, ".bg"))
  bg2 <- import.bedGraph(paste0(sample_name, ".bg"), format="bedGraph")
  bg2 <- sortSeqlevels(bg2)
  bg2 <- sort(bg2)
  seqlengths(bg2) <- seqlengths(txdb)[!grepl("random|chrUn|X|Y|M", names(seqlengths(txdb)))]
  findOverlaps(bg2, reduce(bg2))
  reduce_bg2 <- reduce(bg2)
  export.bw(bg2, paste0(sample_name, ".bw"))
}
```

## Load resulting files to UCSC genome browser for viewing.



