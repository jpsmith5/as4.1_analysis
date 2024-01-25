# A485 ATAC-seq

## 1a. Remove bias from ATAC-seq signals

```console
cd $PROCESSED/A485/
mkdir bam_files/
cd bam_files/
find ../* -type f -name '*dedup.bam' -exec cp {} . \;
export GENOME_FASTA="mm10.fa"
nano sample_list.tsv
```

A485_01
A485_02

```console
$ Get read lengths (looks like 150)
samtools view As4.1_DMSO_01_sort_dedup.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort -n | uniq -c
```

Run seqOutBias
```console
cd /sfs/weka$PROCESSED/A485/bam_files/sob_files/
wget  -O mm10.fa http://awspds.refgenie.databio.org/refgenomes.databio.org/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/fasta__default/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1.fa

seqOutBias mm10.fa As4.1_A485_01_sort_dedup.bam --skip-bed \
	--no-scale --bw=$PROCESSED/A485/bam_files/A485_01.bigWig \
    --only-paired --shift-counts --read-size=150
    
seqOutBias mm10.fa As4.1_A485_02_sort_dedup.bam --skip-bed \
	--no-scale --bw=$PROCESSED/A485/bam_files/A485_02.bigWig \
    --only-paired --shift-counts --read-size=150
```

Copy signal tracks to public location for UCSC usage
```console
cd $PROCESSED/A485/bam_files/
cp *.bigWig $WWW/genome_browser/trackHub/mm10
```

## 1b. Merge replicates and call peaks on the merged samples

```console
cd $PROCESSED/A485/bam_files/

samtools merge -@ 40 /scratch/jps3dp/processed/pepatac/gomez/robert/A485/bam_files/As4.1_A485_merged.bam \
 /scratch/jps3dp/processed/pepatac/gomez/robert/A485/bam_files/As4.1_A485_01_sort_dedup.bam \
 /scratch/jps3dp/processed/pepatac/gomez/robert/A485/bam_files/As4.1_A485_02_sort_dedup.bam 

# Use the combined JQ1/H89 DMSO files as control
samtools merge -@ 40 As4.1_DMSO_merged.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/jq1/bam_files/As4.1_DMSO_01_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/jq1/bam_files/As4.1_DMSO_02_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/H89/bam_files/As4.1_DMSO_01_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/H89/bam_files/As4.1_DMSO_02_sort_dedup.bam

# Merge all for this version of consensus peak calling
samtools merge -@ 40 As4.1_merged.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/jq1/bam_files/As4.1_DMSO_01_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/jq1/bam_files/As4.1_DMSO_02_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/H89/bam_files/As4.1_DMSO_01_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/H89/bam_files/As4.1_DMSO_02_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/A485/bam_files/As4.1_A485_01_sort_dedup.bam \
/scratch/jps3dp/processed/pepatac/gomez/robert/A485/bam_files/As4.1_A485_02_sort_dedup.bam 

macs3 callpeak -t As4.1_A485_merged.bam -f BAMPE -n A485_merged --outdir A485_macs \
	-g mm -B --call-summits --keep-dup 50 -q 0.05 -m 10 200

macs3 callpeak -t As4.1_DMSO_merged.bam -f BAMPE -n DMSO_merged --outdir DMSO_macs \
	-g mm -B --call-summits --keep-dup 50 -q 0.05 -m 10 200

macs3 callpeak -t As4.1_merged.bam -f BAMPE -n As4.1_merged --outdir MERGED_macs \
	-g mm -B --call-summits --keep-dup 50 -q 0.05 -m 10 200

cd A485_macs/
awk '{OFS="\t";} {print $1,$2-100,$3+100,$4,$5}' A485_merged_summits.bed > ../A485_merged_summit_window.bed

cd ../DMSO_macs/
awk '{OFS="\t";} {print $1,$2-100,$3+100,$4,$5}' DMSO_merged_summits.bed > ../DMSO_merged_summit_window.bed
cd ../

cd MERGED_macs/
awk '{OFS="\t";} {print $1,$2-100,$3+100,$4,$5}' As4.1_merged_summits.bed > ../As4.1_merged_summit_window.bed
cd ../

#Example of only using peaks with small p values to keep each set under 200,000
awk '{ if ($5 > 6) { print } }' A485_merged_summit_window.bed > As4.1_A485_summit_window_pval_smallpval.bed

awk '{ if ($5 > 6) { print } }' DMSO_merged_summit_window.bed > As4.1_DMSO_summit_window_pval_smallpval.bed

awk '{ if ($5 > 6) { print } }' As4.1_merged_summit_window.bed > As4.1_summit_window_pval_smallpval.bed
```

## 1c. Merging and removing blacklisted sites

```console
#Merging and removing blacklisted sites
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
gunzip mm10.blacklist.bed.gz
#You can do this with the original files or the small pval files
cat As4.1_A485_summit_window_pval_smallpval.bed | sort -k1,1 -k2,2n | \
	awk ' $2 >= 0 ' | \
	mergeBed -i stdin | \
	intersectBed -wa -v -a stdin -b mm10.blacklist.bed > As4.1_A485_smallpval_rmBlacklist.bed

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
samtools index -@ 8 As4.1_A485_merged.bam
bamCoverage -b As4.1_A485_merged.bam --normalizeUsing BPM --ignoreDuplicates --centerReads --numberOfProcessors max/2 --blackListFileName mm10_blacklist_nonOverlapping.bed --outFileName As4.1_A485_merged.bw

samtools index -@ 8 As4.1_DMSO_merged.bam
bamCoverage -b As4.1_DMSO_merged.bam --normalizeUsing BPM --ignoreDuplicates --centerReads --numberOfProcessors max/2 --blackListFileName mm10_blacklist_nonOverlapping.bed --outFileName As4.1_DMSO_merged.bw
```