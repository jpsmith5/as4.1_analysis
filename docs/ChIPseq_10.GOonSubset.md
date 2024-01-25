### 10. Test for overrepresented GO/KEGG/reactome terms/paths in the "true" set

#### 10a. H3K27Ac
```R
project_name <- "H3K27Ac"
```

##### H89 enriched regions

```R
# GO
true_peaks_H3K27Ac_H89_GO <- ChIPpeakAnno::getEnrichedGO(
    true_peaks_H3K27Ac_H89,
    orgAnn="org.Mm.eg.db",
    maxP=.01,
    minGOterm=10,
    multiAdjMethod="BH",
    condense=TRUE)

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_H89_ChIPpeakAnno_GO.svg"))
enrichmentPlot(true_peaks_H3K27Ac_H89_GO)
dev.off()

# Write GO results to file
for (go in names(true_peaks_H3K27Ac_H89_GO)) {
    file_name <- paste0(project_name, "_DB/", project_name,
                        "_true_peaks_H89_", toupper(go), ".csv")
    setorder(true_peaks_H3K27Ac_H89_GO[[go]], cols = "BH.adjusted.p.value")   
    fwrite(true_peaks_H3K27Ac_H89_GO[[go]], file_name, sep=",")
}

# KEGG
true_peaks_H3K27Ac_H89_KEGG <- ChIPpeakAnno::getEnrichedPATH(
    true_peaks_H3K27Ac_H89,
    "org.Mm.eg.db",
    "KEGGREST",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_H89_ChIPpeakAnno_KEGG.svg"))
enrichmentPlot(true_peaks_H3K27Ac_H89_KEGG)
dev.off()

true_peaks_H3K27Ac_H89_KEGG <- as.data.table(true_peaks_H3K27Ac_H89_KEGG)
setorder(true_peaks_H3K27Ac_H89_KEGG, cols = "BH.adjusted.p.value")  

# KEGG findings
fwrite(true_peaks_H3K27Ac_H89_KEGG,
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_H89_KEGG.csv"))

# REACTOME pathways
true_peaks_H3K27Ac_H89_reactome <- ChIPpeakAnno::getEnrichedPATH(
    true_peaks_H3K27Ac_H89,
    "org.Mm.eg.db",
    "reactome.db",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_H89_ChIPpeakAnno_reactome.svg"))
enrichmentPlot(true_peaks_H3K27Ac_H89_reactome)
dev.off()

true_peaks_H3K27Ac_H89_reactome <- as.data.table(true_peaks_H3K27Ac_H89_reactome)
setorder(true_peaks_H3K27Ac_H89_reactome, cols = "BH.adjusted.p.value")  

# REACTOME pathways: preview data
fwrite(true_peaks_H3K27Ac_H89_reactome,
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_H89_reactome.csv"))
```

##### DMSO enriched regions

```R
# GO
true_peaks_H3K27Ac_DMSO_GO <- ChIPpeakAnno::getEnrichedGO(
    true_peaks_H3K27Ac_DMSO,
    orgAnn="org.Mm.eg.db",
    maxP=.01,
    minGOterm=10,
    multiAdjMethod="BH",
    condense=TRUE)

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_DMSO_ChIPpeakAnno_GO.svg"))
enrichmentPlot(true_peaks_H3K27Ac_DMSO_GO)
dev.off()

# Write GO results to file
for (go in names(true_peaks_H3K27Ac_DMSO_GO)) {
    file_name <- paste0(project_name, "_DB/", project_name,
                        "_true_peaks_DMSO_", toupper(go), ".csv")
    setorder(true_peaks_H3K27Ac_DMSO_GO[[go]], cols = "BH.adjusted.p.value")   
    fwrite(true_peaks_H3K27Ac_DMSO_GO[[go]], file_name, sep=",")
}

# KEGG
true_peaks_H3K27Ac_DMSO_KEGG <- ChIPpeakAnno::getEnrichedPATH(
    true_peaks_H3K27Ac_DMSO,
    "org.Mm.eg.db",
    "KEGGREST",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_DMSO_ChIPpeakAnno_KEGG.svg"))
enrichmentPlot(true_peaks_H3K27Ac_DMSO_KEGG)
dev.off()

true_peaks_H3K27Ac_DMSO_KEGG <- as.data.table(true_peaks_H3K27Ac_DMSO_KEGG)
setorder(true_peaks_H3K27Ac_DMSO_KEGG, cols = "BH.adjusted.p.value")  

# KEGG findings
fwrite(true_peaks_H3K27Ac_DMSO_KEGG,
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_DMSO_KEGG.csv"))

# REACTOME pathways
true_peaks_H3K27Ac_DMSO_reactome <- ChIPpeakAnno::getEnrichedPATH(
    true_peaks_H3K27Ac_DMSO,
    "org.Mm.eg.db",
    "reactome.db",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_DMSO_ChIPpeakAnno_reactome.svg"))
enrichmentPlot(true_peaks_H3K27Ac_DMSO_reactome)
dev.off()

true_peaks_H3K27Ac_DMSO_reactome <- as.data.table(true_peaks_H3K27Ac_DMSO_reactome)
setorder(true_peaks_H3K27Ac_DMSO_reactome, cols = "BH.adjusted.p.value")  

# REACTOME pathways: preview data
fwrite(true_peaks_H3K27Ac_DMSO_reactome,
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_DMSO_reactome.csv"))
```

#### 10b. P300
```R
project_name <- "P300"
```

##### H89 enriched regions

```R
# GO
true_peaks_P300_H89_GO <- ChIPpeakAnno::getEnrichedGO(
    true_peaks_P300_H89,
    orgAnn="org.Mm.eg.db",
    maxP=.01,
    minGOterm=10,
    multiAdjMethod="BH",
    condense=TRUE)

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_H89_ChIPpeakAnno_GO.svg"))
enrichmentPlot(true_peaks_P300_H89_GO)
dev.off()

# Write GO results to file
for (go in names(true_peaks_P300_H89_GO)) {
    file_name <- paste0(project_name, "_DB/", project_name,
                        "_true_peaks_H89_", toupper(go), ".csv")
    setorder(true_peaks_P300_H89_GO[[go]], cols = "BH.adjusted.p.value")   
    fwrite(true_peaks_P300_H89_GO[[go]], file_name, sep=",")
}

# KEGG
true_peaks_P300_H89_KEGG <- ChIPpeakAnno::getEnrichedPATH(
    true_peaks_P300_H89,
    "org.Mm.eg.db",
    "KEGGREST",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_H89_ChIPpeakAnno_KEGG.svg"))
enrichmentPlot(true_peaks_P300_H89_KEGG)
dev.off()

true_peaks_P300_H89_KEGG <- as.data.table(true_peaks_P300_H89_KEGG)
setorder(true_peaks_P300_H89_KEGG, cols = "BH.adjusted.p.value")  

# KEGG findings
fwrite(true_peaks_P300_H89_KEGG,
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_H89_KEGG.csv"))

# REACTOME pathways
true_peaks_P300_H89_reactome <- ChIPpeakAnno::getEnrichedPATH(
    true_peaks_P300_H89,
    "org.Mm.eg.db",
    "reactome.db",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_H89_ChIPpeakAnno_reactome.svg"))
enrichmentPlot(true_peaks_P300_H89_reactome)
dev.off()

true_peaks_P300_H89_reactome <- as.data.table(true_peaks_P300_H89_reactome)
setorder(true_peaks_P300_H89_reactome, cols = "BH.adjusted.p.value")  

# REACTOME pathways: preview data
fwrite(true_peaks_P300_H89_reactome,
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_H89_reactome.csv"))
```

##### DMSO enriched regions

```R
# GO
true_peaks_P300_DMSO_GO <- ChIPpeakAnno::getEnrichedGO(
    true_peaks_P300_DMSO,
    orgAnn="org.Mm.eg.db",
    maxP=.01,
    minGOterm=10,
    multiAdjMethod="BH",
    condense=TRUE)

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_DMSO_ChIPpeakAnno_GO.svg"))
enrichmentPlot(true_peaks_P300_DMSO_GO)
dev.off()

# Write GO results to file
for (go in names(true_peaks_P300_DMSO_GO)) {
    file_name <- paste0(project_name, "_DB/", project_name,
                        "_true_peaks_DMSO_", toupper(go), ".csv")
    setorder(true_peaks_P300_DMSO_GO[[go]], cols = "BH.adjusted.p.value")   
    fwrite(true_peaks_P300_DMSO_GO[[go]], file_name, sep=",")
}

# KEGG
true_peaks_P300_DMSO_KEGG <- ChIPpeakAnno::getEnrichedPATH(
    true_peaks_P300_DMSO,
    "org.Mm.eg.db",
    "KEGGREST",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_DMSO_ChIPpeakAnno_KEGG.svg"))
enrichmentPlot(true_peaks_P300_DMSO_KEGG)
dev.off()

true_peaks_P300_DMSO_KEGG <- as.data.table(true_peaks_P300_DMSO_KEGG)
setorder(true_peaks_P300_DMSO_KEGG, cols = "BH.adjusted.p.value")  

# KEGG findings
fwrite(true_peaks_P300_DMSO_KEGG,
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_DMSO_KEGG.csv"))

# REACTOME pathways
true_peaks_P300_DMSO_reactome <- ChIPpeakAnno::getEnrichedPATH(
    true_peaks_P300_DMSO,
    "org.Mm.eg.db",
    "reactome.db",
    maxP=.01,
    multiAdjMethod="BH")

svglite(paste0(project_name, "_DB/", project_name,
               "_true_peaks_DMSO_ChIPpeakAnno_reactome.svg"))
enrichmentPlot(true_peaks_P300_DMSO_reactome)
dev.off()

true_peaks_P300_DMSO_reactome <- as.data.table(true_peaks_P300_DMSO_reactome)
setorder(true_peaks_P300_DMSO_reactome, cols = "BH.adjusted.p.value")  

# REACTOME pathways: preview data
fwrite(true_peaks_P300_DMSO_reactome,
       paste0(project_name, "_DB/", project_name,
              "_true_peaks_DMSO_reactome.csv"))
```