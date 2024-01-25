### 7. Functional analysis v1 (identifying overrepresented GO/KEGG/reactome terms)

#### 7a. H3K27Ac

```R
project_name <- "H3K27Ac"

# Peak distribution over genomic features
dba_H3K27Ac_peaks_featuresDist <- ChIPpeakAnno::assignChromosomeRegion(
    dba_H3K27Ac_peaksAnno,
    nucleotideLevel=FALSE,
    precedence=c("Promoters", "immediateDownstream", "fiveUTRs",
                 "threeUTRs","Exons", "Introns"),
    TxDb=txdb)

svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_features_distribution.svg"))
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(dba_H3K27Ac_peaks_featuresDist$percentage, las=1, horiz=T)
dev.off()

# Test for overrepresented GO Terms
dba_H3K27Ac_peaks_GO <- ChIPpeakAnno::getEnrichedGO(
    dba_H3K27Ac_peaksAnno,
    orgAnn="org.Mm.eg.db",
    maxP=.01,
    minGOterm=10,
    multiAdjMethod="BH",
    condense=TRUE)

svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_ChIPpeakAnno_GO.svg"))
enrichmentPlot(dba_H3K27Ac_peaks_GO)
dev.off()

# Write GO results to file
for (go in names(dba_H3K27Ac_peaks_GO)) {
    file_name <- paste0(project_name, "_DB/", project_name,
                        "_H89-v-DMSO_", toupper(go), ".csv")
    setorder(dba_H3K27Ac_peaks_GO[[go]], cols = "BH.adjusted.p.value")   
    fwrite(dba_H3K27Ac_peaks_GO[[go]], file_name, sep=",")
}

# KEGG
dba_H3K27Ac_peaks_KEGG <- ChIPpeakAnno::getEnrichedPATH(
    dba_H3K27Ac_peaksAnno,
    "org.Mm.eg.db",
    "KEGGREST",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_ChIPpeakAnno_KEGG.svg"))
enrichmentPlot(dba_H3K27Ac_peaks_KEGG)
dev.off()

dba_H3K27Ac_peaks_KEGG <- as.data.table(dba_H3K27Ac_peaks_KEGG)
setorder(dba_H3K27Ac_peaks_KEGG, cols = "BH.adjusted.p.value")  

# KEGG findings
fwrite(dba_H3K27Ac_peaks_KEGG,
       paste0(project_name, "_DB/", project_name,
              "_H89-v-DMSO_KEGG.csv"))

# REACTOME pathways
dba_H3K27Ac_peaks_pathways <- ChIPpeakAnno::getEnrichedPATH(
    dba_H3K27Ac_peaksAnno,
    "org.Mm.eg.db",
    "reactome.db",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_ChIPpeakAnno_reactome.svg"))
enrichmentPlot(dba_H3K27Ac_peaks_pathways)
dev.off()

dba_H3K27Ac_peaks_pathways <- as.data.table(dba_H3K27Ac_peaks_pathways)
setorder(dba_H3K27Ac_peaks_pathways, cols = "BH.adjusted.p.value")  

# REACTOME pathways: preview data
fwrite(dba_H3K27Ac_peaks_pathways,
       paste0(project_name, "_DB/", project_name,
              "_H89-v-DMSO_reactome.csv"))
```

#### 7b. P300

```R
project_name <- "P300"

# Peak distribution over genomic features
dba_P300_peaks_featuresDist <- ChIPpeakAnno::assignChromosomeRegion(
    dba_P300_peaksAnno,
    nucleotideLevel=FALSE,
    precedence=c("Promoters", "immediateDownstream", "fiveUTRs",
                 "threeUTRs","Exons", "Introns"),
    TxDb=txdb)

svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_features_distribution.svg"))
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(dba_P300_peaks_featuresDist$percentage, las=1, horiz=T)
dev.off()

# Test for overrepresented GO Terms
dba_P300_peaks_GO <- ChIPpeakAnno::getEnrichedGO(
    dba_P300_peaksAnno,
    orgAnn="org.Mm.eg.db",
    maxP=.01,
    minGOterm=10,
    multiAdjMethod="BH",
    condense=TRUE)

svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_ChIPpeakAnno_GO.svg"))
enrichmentPlot(dba_P300_peaks_GO)
dev.off()

# Write GO results to file
for (go in names(dba_P300_peaks_GO)) {
    file_name <- paste0(project_name, "_DB/", project_name,
                        "_H89-v-DMSO_", toupper(go), ".csv")
    setorder(dba_P300_peaks_GO[[go]], cols = "BH.adjusted.p.value")   
    fwrite(dba_P300_peaks_GO[[go]], file_name, sep=",")
}

# KEGG
dba_P300_peaks_KEGG <- ChIPpeakAnno::getEnrichedPATH(
    dba_P300_peaksAnno,
    "org.Mm.eg.db",
    "KEGGREST",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_ChIPpeakAnno_KEGG.svg"))
enrichmentPlot(dba_P300_peaks_KEGG)
dev.off()

dba_P300_peaks_KEGG <- as.data.table(dba_P300_peaks_KEGG)
setorder(dba_P300_peaks_KEGG, cols = "BH.adjusted.p.value")  

# KEGG findings
fwrite(dba_P300_peaks_KEGG,
       paste0(project_name, "_DB/", project_name,
              "_H89-v-DMSO_KEGG.csv"))

# REACTOME pathways
dba_P300_peaks_pathways <- ChIPpeakAnno::getEnrichedPATH(
    dba_P300_peaksAnno,
    "org.Mm.eg.db",
    "reactome.db",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_ChIPpeakAnno_reactome.svg"))
enrichmentPlot(dba_P300_peaks_pathways)
dev.off()

dba_P300_peaks_pathways <- as.data.table(dba_P300_peaks_pathways)
setorder(dba_P300_peaks_pathways, cols = "BH.adjusted.p.value")  

# REACTOME pathways: preview data
fwrite(dba_P300_peaks_pathways,
       paste0(project_name, "_DB/", project_name,
              "_H89-v-DMSO_reactome.csv"))
```