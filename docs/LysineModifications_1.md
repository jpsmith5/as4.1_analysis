# Lysine Modifications

 - H3K27ac and H2BK5ac (both deposited by CBP/p300)
  - performed immunofluorescence (IF) for H3K27ac and H2BK5ac in A485 treated cells. As4.1 cells treated with A485 for 1 hour had a large decrease in H3K27ac and H2BK5ac intensity by IF- the marks were undetectable compared to DMSO treated controls. 

## 1a. Set up sample_table.csv

`sample_table.csv`:
```
group	replicate	fastq_1	fastq_2	control
H3K4me1	1	H3K4me1_1_CKDL210025996-1a_H7VKCDSX3_L2_1.fq.gz	H3K4me1_1_CKDL210025996-1a_H7VKCDSX3_L2_2.fq.gz	igg_ctrl
H3K4me1	2	H3K4me1_2_CKDL210026002-1a_H3H55DSX3_L2_1.fq.gz	H3K4me1_2_CKDL210026002-1a_H3H55DSX3_L2_2.fq.gz	igg_ctrl
H3K4me3	1	H3K4me3_1_CKDL210025997-1a_H7VKCDSX3_L2_1.fq.gz	H3K4me3_1_CKDL210025997-1a_H7VKCDSX3_L2_2.fq.gz	igg_ctrl
H3K4me3	2	H3K4me3_2_CKDL210026003-1a_H3H55DSX3_L2_1.fq.gz	H3K4me3_2_CKDL210026003-1a_H3H55DSX3_L2_2.fq.gz	igg_ctrl
H2B5Kac	1	H2B5Kac_1_CKDL210025992-1a_H3H55DSX3_L1_1.fq.gz	H2B5Kac_1_CKDL210025992-1a_H3H55DSX3_L1_2.fq.gz	igg_ctrl
H2B5Kac	2	H2B5Kac_2_CKDL210025998-1a_H3H55DSX3_L3_1.fq.gz	H2B5Kac_2_CKDL210025998-1a_H3H55DSX3_L3_2.fq.gz	igg_ctrl
H4K16ac	1	H4K16ac_1_CKDL210025995-1a_H7VKCDSX3_L2_1.fq.gz	H4K16ac_1_CKDL210025995-1a_H7VKCDSX3_L2_2.fq.gz	igg_ctrl
H4K16ac	2	H4K16ac_2_CKDL210026001-1a_H3H55DSX3_L2_1.fq.gz	H4K16ac_2_CKDL210026001-1a_H3H55DSX3_L2_2.fq.gz	igg_ctrl
H4K5ac	1	H4K5ac_1_CKDL210025994-1a_H3H55DSX3_L1_1.fq.gz	H4K5ac_1_CKDL210025994-1a_H3H55DSX3_L1_2.fq.gz	igg_ctrl
H4K5ac	2	H4K5ac_2_CKDL210026000-1a_H3H55DSX3_L3_1.fq.gz	H4K5ac_2_CKDL210026000-1a_H3H55DSX3_L3_2.fq.gz	igg_ctrl
H3K27ac	1	H3K27ac_1_CKDL210025993-1a_H3H55DSX3_L1_1.fq.gz	H3K27ac_1_CKDL210025993-1a_H3H55DSX3_L1_2.fq.gz	igg_ctrl
H3K27ac	2	H3K27ac_2_CKDL210025999-1a_H3H55DSX3_L3_1.fq.gz	H3K27ac_2_CKDL210025999-1a_H3H55DSX3_L3_2.fq.gz	igg_ctrl
igg_ctrl	1	IgG_CKDL210026004-1a_H3H55DSX3_L2_1.fq.gz	IgG_CKDL210026004-1a_H3H55DSX3_L2_2.fq.gz	
```

## 1b. Run pipeline to generate signal

```console
nextflow run nf-core/cutandrun --input $DATA/lysine_modifications/sample_table.csv --outdir $PROCESSED/lysine_modifications --genome GRCm38 -profile singularity -c nextflow_rivanna.conf --macs_gsize 2406655830 -bg -resume --max_cpus=4
--dt_calc_all_matrix FALSE
```

## 1c. Convert consensus peaks to bigBed for UCSC genome browser

```console
cd $PROCESSED/lysine_modifications/consensus_peaks/

awk -F'\t' '{print "chr"$1, $2, $3}' H3K4me1.consensus.peak_counts.bed | sort -k1,1 -k2,2n > H3K4me1_consensus_peaks.bed3
awk -F'\t' '{print "chr"$1, $2, $3}' H3K4me3.consensus.peak_counts.bed | sort -k1,1 -k2,2n  > H3K4me3_consensus_peaks.bed3
awk -F'\t' '{print "chr"$1, $2, $3}' H2B5Kac.consensus.peak_counts.bed | sort -k1,1 -k2,2n  > H2B5Kac_consensus_peaks.bed3
awk -F'\t' '{print "chr"$1, $2, $3}' H3K27ac.consensus.peak_counts.bed | sort -k1,1 -k2,2n  > H3K27ac_consensus_peaks.bed3
awk -F'\t' '{print "chr"$1, $2, $3}' H4K16ac.consensus.peak_counts.bed | sort -k1,1 -k2,2n > H4K16ac_consensus_peaks.bed3
awk -F'\t' '{print "chr"$1, $2, $3}' H4K5ac.consensus.peak_counts.bed | sort -k1,1 -k2,2n > H4K5ac_consensus_peaks.bed3

sed -i 's/chrMT/chrM/g' *.bed3

bedToBigBed -type=bed3 H3K4me1_consensus_peaks.bed3 mm10.chrom.sizes H3K4me1_consensus_peaks.bb
bedToBigBed -type=bed3 H3K4me3_consensus_peaks.bed3 mm10.chrom.sizes H3K4me3_consensus_peaks.bb
bedToBigBed -type=bed3 H2B5Kac_consensus_peaks.bed3 mm10.chrom.sizes H2B5Kac_consensus_peaks.bb
bedToBigBed -type=bed3 H3K27ac_consensus_peaks.bed3 mm10.chrom.sizes H3K27ac_consensus_peaks.bb
bedToBigBed -type=bed3 H4K16ac_consensus_peaks.bed3 mm10.chrom.sizes H4K16ac_consensus_peaks.bb
bedToBigBed -type=bed3 H4K5ac_consensus_peaks.bed3 mm10.chrom.sizes H4K5ac_consensus_peaks.bb

cp *.bb $WWW/genome_browser/trackHub/mm10/
```

## 1d. Add signal tracks

```console
cd $PROCESSED/lysine_modifications/

cp $PROCESSED//lysine_modifications/03_peak_calling/03_bed_to_bigwig/*.bigWig $WWW/genome_browser/trackHub/mm10/

mkdir bigwig_files/
cp $PROCESSED/lysine_modifications/03_peak_calling/03_bed_to_bigwig/*.bigWig $PROCESSED/lysine_modifications/bigwig_files

bigWigMerge H2B5Kac_R1.bigWig H2B5Kac_R2.bigWig H2B5Kac.bedGraph
awk -F'\t' '{print "chr"$1, $2, $3, $4}' H2B5Kac.bedGraph > tmp && mv tmp H2B5Kac.bedGraph
sed -i 's/chrMT/chrM/g' H2B5Kac.bedGraph
#TODO: return here
LC_COLLATE=C sort -k1,1 -k2,2n H2B5Kac.bedGraph > tmp && mv tmp H2B5Kac.bedGraph
bedGraphToBigWig H2B5Kac.bedGraph mm10.chrom.sizes H2B5Kac.bigWig

bigWigMerge H3K27ac_R1.bigWig H3K27ac_R2.bigWig H3K27ac.bedGraph
awk -F'\t' '{print "chr"$1, $2, $3, $4}' H3K27ac.bedGraph > tmp && mv tmp H3K27ac.bedGraph
sed -i 's/chrMT/chrM/g' H3K27ac.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n H3K27ac.bedGraph > tmp && mv tmp H3K27ac.bedGraph
bedGraphToBigWig H3K27ac.bedGraph mm10.chrom.sizes H3K27ac.bigWig

bigWigMerge H3K4me1_R1.bigWig H3K4me1_R2.bigWig H3K4me1.bedGraph
awk -F'\t' '{print "chr"$1, $2, $3, $4}' H3K4me1.bedGraph > tmp && mv tmp H3K4me1.bedGraph
sed -i 's/chrMT/chrM/g' H3K4me1.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n H3K4me1.bedGraph > tmp && mv tmp H3K4me1.bedGraph
bedGraphToBigWig H3K4me1.bedGraph mm10.chrom.sizes H3K4me1.bigWig

bigWigMerge H3K4me3_R1.bigWig H3K4me3_R2.bigWig H3K4me3.bedGraph
awk -F'\t' '{print "chr"$1, $2, $3, $4}' H3K4me3.bedGraph > tmp && mv tmp H3K4me3.bedGraph
sed -i 's/chrMT/chrM/g' H3K4me3.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n H3K4me3.bedGraph > tmp && mv tmp H3K4me3.bedGraph
bedGraphToBigWig H3K4me3.bedGraph mm10.chrom.sizes H3K4me3.bigWig

bigWigMerge H4K16ac_R1.bigWig H4K16ac_R2.bigWig H4K16ac.bedGraph
awk -F'\t' '{print "chr"$1, $2, $3, $4}' H4K16ac.bedGraph > tmp && mv tmp H4K16ac.bedGraph
sed -i 's/chrMT/chrM/g' H4K16ac.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n H4K16ac.bedGraph > tmp && mv tmp H4K16ac.bedGraph
bedGraphToBigWig H4K16ac.bedGraph mm10.chrom.sizes H4K16ac.bigWig

bigWigMerge H4K5ac_R1.bigWig H4K5ac_R2.bigWig H4K5ac.bedGraph
awk -F'\t' '{print "chr"$1, $2, $3, $4}' H4K5ac.bedGraph > tmp && mv tmp H4K5ac.bedGraph
sed -i 's/chrMT/chrM/g' H4K5ac.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n H4K5ac.bedGraph > tmp && mv tmp H4K5ac.bedGraph
bedGraphToBigWig H4K5ac.bedGraph mm10.chrom.sizes H4K5ac.bigWig

cp $PROCESSED/lysine_modifications/bigwig_files/*.bigWig $WWW/genome_browser/trackHub/mm10/
```

