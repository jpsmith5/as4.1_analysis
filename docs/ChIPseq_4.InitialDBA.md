### 4a. H3K27Ac ChIP-seq with H89 or DMSO treated As4.1 cells

```R
dba_H3K27Ac <- calcDBAobject(samples=samples_H3K27Ac, project_name="H3K27Ac")
plotDBAobject(dba_H3K27Ac, project_name="H3K27Ac")
svglite(paste0("H3K27Ac", "_DB/", "H3K27Ac",
               "_H89-v-DMSO_MAplot.svg"))
dba.plotMA(dba_H3K27Ac)
dev.off()

dba_H3K27Ac_DB <- dba.report(dba_H3K27Ac)

# Add gene annotations and look at fold changes of peaks 
# nearest to genes of interest.
## Loading TSS Annotation Obtained From BiomaRt
data(TSS.mouse.GRCm38)

# Annotate peaks with information on closest TSS using precompiled annotation data
dba_H3K27Ac_peaksAnno <- annotatePeakInBatch(dba_H3K27Ac_DB,
                                             AnnotationData=TSS.mouse.GRCm38)

svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_ChIPpeakAnno_PIE.svg"))
pie1(table(as.data.frame(dba_H3K27Ac_peaksAnno)$insideFeature))
dev.off()

dba_H3K27Ac_peaksAnno <- addGeneIDs(dba_H3K27Ac_peaksAnno,
                                    orgAnn="org.Mm.eg.db",
                                    IDs2Add=c("symbol"))
```

### 4b. P300 ChIP-seq with H89 or DMSO treated As4.1 cells

```R
dba_P300 <- calcDBAobject(samples=samples_P300, project_name="P300")
plotDBAobject(dba_P300, project_name="P300")
svglite(paste0("P300", "_DB/", "P300",
               "_H89-v-DMSO_MAplot.svg"))
dba.plotMA(dba_P300)
dev.off()

dba_P300_DB <- dba.report(dba_P300)

# Add gene annotations and look at fold changes of peaks 
# nearest to genes of interest.
## Loading TSS Annotation Obtained From BiomaRt
data(TSS.mouse.GRCm38)

# Annotate peaks with information on closest TSS using precompiled annotation data
dba_P300_peaksAnno <- annotatePeakInBatch(dba_P300_DB,
                                          AnnotationData=TSS.mouse.GRCm38)

svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_ChIPpeakAnno_PIE.svg"))
pie1(table(as.data.frame(dba_P300_peaksAnno)$insideFeature))
dev.off()

dba_P300_peaksAnno <- addGeneIDs(dba_P300_peaksAnno,
                                    orgAnn="org.Mm.eg.db",
                                    IDs2Add=c("symbol"))
dba_P300_peaksAnno[grep("Ren1", dba_P300_peaksAnno$symbol),]
# Ren1 P300 drops with H89 treatment
dba_P300_peaksAnno[grep("Nfix", dba_P300_peaksAnno$symbol),]
# N/A
```