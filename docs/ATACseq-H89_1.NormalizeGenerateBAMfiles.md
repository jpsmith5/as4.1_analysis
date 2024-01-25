# H89 ATAC-seq 

## 1a. Remove bias from ATAC-seq signals

Confirm read lengths

```console
samtools view As4.1_DMSO_03_sort_dedup.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort -n | uniq -c
nano sample_list.tsv
```

DMSO_01
DMSO_02
DMSO_03
H89_01
H89_02

```console
wget -O mm10.fa.tgz http://refgenomes.databio.org/v3/assets/archive/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/fasta?tag=default
tar xvfz mm10.fa.tgz 
export GENOME_FASTA="default/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1.fa"

while read sample
do
	echo $sample
	seqOutBias ${GENOME_FASTA} As4.1_${sample}_sort_dedup.bam --skip-bed \
	--no-scale --bw=${sample}.bigWig --only-paired --shift-counts --read-size=150
done < sample_list.tsv
```

## 1b. Merge replicates and call peaks on the merged samples

```console
samtools merge -@ 24 As4.1_H89_merged.bam \
 As4.1_H89_01_sort_dedup.bam \
 As4.1_H89_02_sort_dedup.bam 
samtools merge -@ 24 As4.1_DMSO_merged.bam \
 As4.1_DMSO_01_sort_dedup.bam \
 As4.1_DMSO_02_sort_dedup.bam \
 As4.1_DMSO_03_sort_dedup.bam
samtools merge -@ 24 As4.1_merged.bam \
 As4.1_DMSO_01_sort_dedup.bam \
 As4.1_DMSO_02_sort_dedup.bam \
 As4.1_DMSO_03_sort_dedup.bam \
 As4.1_H89_01_sort_dedup.bam \
 As4.1_H89_02_sort_dedup.bam 

macs3 callpeak -t As4.1_H89_merged.bam -f BAMPE -n H89_merged --outdir H89_macs \
 -g mm -B --call-summits --keep-dup 50 -q 0.05 -m 10 200

macs3 callpeak -t As4.1_DMSO_merged.bam -f BAMPE -n DMSO_merged --outdir DMSO_macs \
 -g mm -B --call-summits --keep-dup 50 -q 0.05 -m 10 200

macs3 callpeak -t As4.1_merged.bam -f BAMPE -n As4.1_merged --outdir As4.1_macs \
 -g mm -B --call-summits --keep-dup 50 -q 0.05 -m 10 200

cd H89_macs/
awk '{OFS="\t";} {print $1,$2-100,$3+100,$4,$5}' H89_merged_summits.bed > ../H89_merged_summit_window.bed

cd ../DMSO_macs/
awk '{OFS="\t";} {print $1,$2-100,$3+100,$4,$5}' DMSO_merged_summits.bed > ../DMSO_merged_summit_window.bed
cd ../

cd As4.1_macs/
awk '{OFS="\t";} {print $1,$2-100,$3+100,$4,$5}' As4.1_merged_summits.bed > ../As4.1_merged_summit_window.bed
cd ../

# Using peaks with small p values to keep each set at reasonable number (sub 200k)
awk '{ if ($5 > 6) { print } }' H89_merged_summit_window.bed > As4.1_H89_summit_window_pval_smallpval.bed
awk '{ if ($5 > 6) { print } }' DMSO_merged_summit_window.bed > As4.1_DMSO_summit_window_pval_smallpval.bed
awk '{ if ($5 > 6) { print } }' As4.1_merged_summit_window.bed > As4.1_summit_window_pval_smallpval.bed
```

## 1c. Merging and removing blacklisted sites

```console
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
gunzip mm10.blacklist.bed.gz

cat As4.1_H89_summit_window_pval_smallpval.bed | sort -k1,1 -k2,2n | \
	awk ' $2 >= 0 ' | \
	mergeBed -i stdin | \
	intersectBed -wa -v -a stdin -b mm10.blacklist.bed > As4.1_H89_smallpval_rmBlacklist.bed

cat As4.1_DMSO_summit_window_pval_smallpval.bed | sort -k1,1 -k2,2n | \
	awk ' $2 >= 0 ' | \
	mergeBed -i stdin | \
	intersectBed -wa -v -a stdin -b mm10.blacklist.bed > As4.1_DMSO_smallpval_rmBlacklist.bed
    
cat As4.1_summit_window_pval_smallpval.bed | sort -k1,1 -k2,2n | \
	awk ' $2 >= 0 ' | \
	mergeBed -i stdin | \
	intersectBed -wa -v -a stdin -b mm10.blacklist.bed > As4.1_smallpval_rmBlacklist.bed
```

## 1d. Convert BAM to bigWig for visualization
 - must reduce the minimum interval between blacklisted regions to satisfy bamCoverage
 
```R
library(data.table)
library(GenomicRanges)
library(rtracklayer)

blacklist <- fread('mm10.blacklist.bed')
colnames(blacklist) <- c("chr", "start", "end")
blacklist_GR <- makeGRangesFromDataFrame(blacklist)
blacklist_GR_final <- reduce(blacklist_GR, min.gapwidth=160L)
export.bed(blacklist_GR_final, 'mm10_blacklist_nonOverlapping.bed')
```

```console
samtools index -@ 8 As4.1_H89_merged.bam
bamCoverage -b As4.1_H89_merged.bam --normalizeUsing BPM --ignoreDuplicates \
 --centerReads --numberOfProcessors max/2 \
 --blackListFileName mm10_blacklist_nonOverlapping.bed \
 --outFileName As4.1_H89_merged.bw

samtools index -@ 8 As4.1_DMSO_merged.bam
bamCoverage -b As4.1_DMSO_merged.bam --normalizeUsing BPM --ignoreDuplicates \
 --centerReads --numberOfProcessors max/2 \
 --blackListFileName mm10_blacklist_nonOverlapping.bed \
 --outFileName As4.1_DMSO_merged.bw
```
