# JQ1 ATAC-seq

## 1a. Remove bias from ATAC-seq signals

```console
cd $PROCESSED/jq1/
mkdir bam_files/
cd bam_files/
find ../* -type f -name '*dedup.bam' -exec cp {} . \;

mv As4.1_Jq1_minus_01_sort_dedup.bam As4.1_DMSO_01_sort_dedup.bam
mv As4.1_Jq1_minus_02_sort_dedup.bam As4.1_DMSO_02_sort_dedup.bam
mv As4.1_Jq1_plus_01_sort_dedup.bam As4.1_JQ1_01_sort_dedup.bam
mv As4.1_Jq1_plus_02_sort_dedup.bam As4.1_JQ1_02_sort_dedup.bam

export GENOME_FASTA="mm10.fa"

nano sample_list.tsv
```

DMSO_01
DMSO_02
JQ1_01
JQ1_02

```console
$ Get read lengths (looks like 150)
samtools view As4.1_DMSO_01_sort_dedup.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort -n | uniq -c
```

Run seqOutBias on cluster

```console
samtools index -@ 40 As4.1_JQ1_01_sort_dedup.bam
samtools index -@ 40 As4.1_JQ1_02_sort_dedup.bam
samtools index -@ 40 As4.1_DMSO_01_sort_dedup.bam
samtools index -@ 40 As4.1_DMSO_02_sort_dedup.bam

while read sample
do
	echo $sample
	seqOutBias ${GENOME_FASTA} As4.1_${sample}_sort_dedup.bam --skip-bed \
	--no-scale --bw=${sample}.bigWig --only-paired --shift-counts --read-size=150
done < sample_list.tsv
```

Copy signal tracks to public location for UCSC usage
```console
cd $PROCESSED/jq1/bam_files/
cp *.bigWig $WWW/genome_browser/trackHub/mm10
```

## 1b. Merge replicates and call peaks on the merged samples

```console
cd $PROCESSED/jq1/bam_files/

samtools merge -@ 40 As4.1_JQ1_merged.bam As4.1_JQ1_01_sort_dedup.bam As4.1_JQ1_02_sort_dedup.bam 
# Use the combined JQ1/H89 DMSO files as control
samtools merge -@ 40 As4.1_DMSO_merged.bam \
$PROCESSED/jq1/bam_files/As4.1_DMSO_01_sort_dedup.bam \
$PROCESSED/jq1/bam_files/As4.1_DMSO_02_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/H89/bam_files/As4.1_DMSO_01_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/H89/bam_files/As4.1_DMSO_02_sort_dedup.bam

# Merge all for this version of consensus peak calling
samtools merge -@ 40 As4.1_merged.bam \
$PROCESSED/jq1/bam_files/As4.1_DMSO_01_sort_dedup.bam \
$PROCESSED/jq1/bam_files/As4.1_DMSO_02_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/H89/bam_files/As4.1_DMSO_01_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/H89/bam_files/As4.1_DMSO_02_sort_dedup.bam \
$PROCESSED/jq1/bam_files/As4.1_JQ1_01_sort_dedup.bam \
$PROCESSED/jq1/bam_files/As4.1_JQ1_02_sort_dedup.bam 

macs3 callpeak -t As4.1_JQ1_merged.bam -f BAMPE -n JQ1_merged --outdir JQ1_macs \
	-g mm -B --call-summits --keep-dup 50 -q 0.05 -m 10 200

macs3 callpeak -t As4.1_DMSO_merged.bam -f BAMPE -n DMSO_merged --outdir DMSO_macs \
	-g mm -B --call-summits --keep-dup 50 -q 0.05 -m 10 200

macs3 callpeak -t As4.1_merged.bam -f BAMPE -n As4.1_merged --outdir As4.1_macs \
	-g mm -B --call-summits --keep-dup 50 -q 0.05 -m 10 200

cd JQ1_macs/
awk '{OFS="\t";} {print $1,$2-100,$3+100,$4,$5}' JQ1_merged_summits.bed > ../JQ1_merged_summit_window.bed

cd ../DMSO_macs/
awk '{OFS="\t";} {print $1,$2-100,$3+100,$4,$5}' DMSO_merged_summits.bed > ../DMSO_merged_summit_window.bed
cd ../

cd As4.1_macs/
awk '{OFS="\t";} {print $1,$2-100,$3+100,$4,$5}' As4.1_merged_summits.bed > ../As4.1_merged_summit_window.bed
cd ../

#Example of only using peaks with small p values to keep each set under 200,000
awk '{ if ($5 > 6) { print } }' JQ1_merged_summit_window.bed > As4.1_JQ1_summit_window_pval_smallpval.bed

awk '{ if ($5 > 6) { print } }' DMSO_merged_summit_window.bed > As4.1_DMSO_summit_window_pval_smallpval.bed

awk '{ if ($5 > 6) { print } }' As4.1_merged_summit_window.bed > As4.1_summit_window_pval_smallpval.bed
```

## 1c. Merging and removing blacklisted sites

```console
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
gunzip mm10.blacklist.bed.gz
#You can do this with the original files or the small pval files
cat As4.1_JQ1_summit_window_pval_smallpval.bed | sort -k1,1 -k2,2n | \
	awk ' $2 >= 0 ' | \
	mergeBed -i stdin | \
	intersectBed -wa -v -a stdin -b mm10.blacklist.bed > As4.1_JQ1_smallpval_rmBlacklist.bed

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
q()
```

```console
samtools index -@ 8 As4.1_JQ1_merged.bam
bamCoverage -b As4.1_JQ1_merged.bam --normalizeUsing BPM --ignoreDuplicates --centerReads --numberOfProcessors max/2 --blackListFileName mm10_blacklist_nonOverlapping.bed --outFileName As4.1_JQ1_merged.bw

samtools index -@ 8 As4.1_DMSO_merged.bam
bamCoverage -b As4.1_DMSO_merged.bam --normalizeUsing BPM --ignoreDuplicates --centerReads --numberOfProcessors max/2 --blackListFileName mm10_blacklist_nonOverlapping.bed --outFileName As4.1_DMSO_merged.bw

samtools index -@ 8 As4.1_merged.bam
bamCoverage -b As4.1_merged.bam --normalizeUsing BPM --ignoreDuplicates --centerReads --numberOfProcessors max/2 --blackListFileName mm10_blacklist_nonOverlapping.bed --outFileName As4.1_merged.bw
```