# ChIP-seq Analysis

 - ChIP-seq for H89(Protein Kinase A) and DMSO treated As4.1 cells
  - H3K27ac
  - p300

## 1. H3K27Ac Preprocessing and intial pipeline analysis

raw data files:
```console
1_08HL_015CUVirginia_DMSO-1_H3K27Ac_mm-dm_i89.fastq.gz  
2_08HM_015CUVirginia_DMSO-2_H3K27Ac_mm-dm_i90.fastq.gz  
3_08HN_015CUVirginia_H89-1_H3K27Ac_mm-dm_i91.fastq.gz
4_08HO_015CUVirginia_H89-2_H3K27Ac_mm-dm_i92.fastq.gz
5_08I8_015CUVirginia_Pooled_Input_mm_i54.fastq.gz
```

 - read_length: 75
 - macs effective genome size: 2406655830

```console
nextflow run nf-core/chipseq --input $DATA/ChIP-seq/metadata/H89_samplesheet.csv --outdir $PROCESSED/chipseq/H89 --genome GRCm38 -profile singularity -c nextflow_rivanna.conf --macs_gsize 2406655830 -bg -resume --max_cpus=1
```

## 2. P300 Preprocessing and intial pipeline analysis

raw data files:
```
1_0BKX_01P2UVirginia_DMSO-1_p300_mm_i87.fastq.gz
2_0BKY_01P2UVirginia_DMSO-2_p300_mm_i88.fastq.gz
3_0BKZ_01P2UVirginia_H89-1_p300_mm_i90.fastq.gz
4_0BL0_01P2UVirginia_H89-2_p300_mm_i91.fastq.gz
5_08I8_015CUVirginia_Pooled_Input_mm_i54.fastq.gz

nextflow run nf-core/chipseq --input $DATA/ChIP-seq/metadata/P300_samplesheet.csv --outdir $PROCESSED/chipseq/P300 --genome GRCm38 -profile singularity -c /project/gomezlab/code/nextflow_rivanna.conf --macs_gsize 2406655830 -bg -resume --max_cpus=1
```