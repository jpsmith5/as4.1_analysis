## 9. Move files to public genome_browser/ location

Format and move scRNA-seq signal tracks to genome_browser/ locations

```console
export PROJECT_WWW="$WWW/genome_browser/trackHub/mm10/"

bamCoverage --numberOfProcessors 40 --effectiveGenomeSize 2620345972 --normalizeUsing RPKM -b $PROCESSED/scRNA/As4-1-H89_scRNAseq/outs/possorted_genome_bam.bam -o ${PROJECT_WWW}/As4.1_H89_scRNAseq_RPKM.bigWig

bamCoverage --numberOfProcessors 40 --effectiveGenomeSize 2620345972 --normalizeUsing RPKM -b $PROCESSED/scRNA/As4-1-control_scRNAseq/outs/possorted_genome_bam.bam -o ${PROJECT_WWW}/As4.1_DMSO_scRNAseq_RPKM.bigWig

```
