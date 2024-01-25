## 5. Perform Gene Set Enrichment Analysis (GSEA) for each population of cells

### 5a. Export expression matrix for each data set prior to GSEA

```R
ctrl_sct <- combined_sct[, combined_sct$orig.ident == "ctrl"]
h89_sct  <- combined_sct[, combined_sct$orig.ident == "h89"]

ctrl_GEX <- GetAssayData(object = ctrl_sct[["SCT"]], slot = "scale.data")
h89_GEX  <- GetAssayData(object = h89_sct[["SCT"]], slot = "scale.data")
```

### 5b. Update palette

Update the color palette with the cluster names for consistency in plotting

```R
pal_named        <- copy(pal)
names(pal_named) <- unname(new_cluster_ids)
pal_named[['remodel']] <- "#C90076"
```

### 5c. Set up ontologies to evaluate

[Download GO files](https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp?targetSpeciesDB=Mouse#MH)

```R
GO_hallmarks     <- "msigdb/mh.all.v2023.1.Mm.symbols.gmt"
GO_canonical     <- "msigdb/m2.all.v2023.1.Mm.symbols.gmt"
GO_TF_targets    <- "msigdb/m3.all.v2023.1.Mm.symbols.gmt"
GO_same_ontology <- "msigdb/m5.all.v2023.1.Mm.symbols.gmt"
GO_cell_sig      <- "msigdb/m8.all.v2023.1.Mm.symbols.gmt"

GO_list <- c(GO_hallmarks, GO_canonical, GO_TF_targets,
             GO_same_ontology, GO_cell_sig)
names(GO_list) <- c("hallmarks", "canonical", "TF_targets", 
                    "same_ontology", "cell_signature")
```

### 5d. Plot enriched ontologies in control versus treated cells

```R
for (GO_file in names(GO_list)) {
    res = GSEA(ctrl_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("As4.1_CTRL_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

for (GO_file in names(GO_list)) {
    res = GSEA(h89_gene_list, GO_list[GO_file], pval = 0.05)
    svglite(paste0("As4.1_H89_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}
```

### 5e. Evaluate individual cell populations identified earlier

```R
named_sct       <- PrepSCTFindMarkers(named_sct, assay="SCT")
cluster_markers <- FindAllMarkers(named_sct, only.pos = FALSE, assay="SCT",
                               min.pct = 0.25, logfc.threshold = 0.25,
                               test.use="MAST", random.seed=99)
```

#### Cluster 0 == all+ 

This cluster had moderate expression of markers of each additional cluster.

```R
cluster0_genes <- cluster_markers[cluster_markers$cluster=="all+",]$avg_log2FC
names(cluster0_genes) = cluster_markers[cluster_markers$cluster=="all+",]$gene
cluster0_genes = sort(cluster0_genes, decreasing = TRUE)
cluster0_genes = cluster0_genes[!duplicated(names(cluster0_genes))]
head(cluster0_genes)

for (GO_file in names(GO_list)) {
    res = GSEA(cluster0_genes, GO_list[GO_file], pval = 0.05)
    svglite(paste0("As4.1_all+_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

cluster0 <- named_sct[, named_sct$named_cluster %in% "all+"] 
table(cluster0$orig.ident)
```

| cluster | ctrl (N) | H89 (N) |
|---------|----------|---------|
| all+    | 1512     | 1756    |


#### Cluster 1 == remodel

```R
cluster1_genes <- cluster_markers[cluster_markers$cluster=="remodel",]$avg_log2FC
names(cluster1_genes) = cluster_markers[cluster_markers$cluster=="remodel",]$gene
cluster1_genes = sort(cluster1_genes, decreasing = TRUE)
cluster1_genes = cluster1_genes[!duplicated(names(cluster1_genes))]
head(cluster1_genes)

for (GO_file in names(GO_list)) {
    res = GSEA(cluster1_genes, GO_list[GO_file], pval = 0.05)
    svglite(paste0("As4.1_remodel_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

cluster1 <- named_sct[, named_sct$named_cluster %in% "remodel"] 
table(cluster1$orig.ident)
```

| cluster | ctrl (N) | H89 (N) |
|---------|----------|---------|
| remodel | 562      | 728     |

 - G2M checkpoint is UP
 - E2F targets UP

Looks like this is a cluster of actively dividing / remodeling cells


#### Cluster 2 == Khdrbs3|Acaa1b+

Marked by the highest and near exclusive expression of *Khdrbs3* and *Acaa1b*.

```R
cluster2_genes <- cluster_markers[cluster_markers$cluster=="Khdrbs3|Acaa1b+",]$avg_log2FC
names(cluster2_genes) = cluster_markers[cluster_markers$cluster=="Khdrbs3|Acaa1b+",]$gene
cluster2_genes = sort(cluster2_genes, decreasing = TRUE)
cluster2_genes = cluster2_genes[!duplicated(names(cluster2_genes))]
head(cluster2_genes)

for (GO_file in names(GO_list)) {
    res = GSEA(cluster2_genes, GO_list[GO_file], pval = 0.05)
    svglite(paste0("As4.1_Khdrbs3-Acaa1b+_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

cluster2 <- named_sct[, named_sct$named_cluster %in% "Khdrbs3|Acaa1b+"] 
table(cluster2$orig.ident)
```

| cluster         | ctrl (N) | H89 (N) |
|-----------------|----------|---------|
| Khdrbs3|Acaa1b+ | 588      | 595     |

EMT is most DOWNregulated GO term.

#### Cluster 3 == Ren1_LO

```R
Ren1_LO_genes <- cluster_markers[cluster_markers$cluster=="Ren1_LO",]$avg_log2FC
names(Ren1_LO_genes) = cluster_markers[cluster_markers$cluster=="Ren1_LO",]$gene
Ren1_LO_genes = sort(Ren1_LO_genes, decreasing = TRUE)
Ren1_LO_genes = Ren1_LO_genes[!duplicated(names(Ren1_LO_genes))]
head(Ren1_LO_genes)

for (GO_file in names(GO_list)) {
    res = GSEA(Ren1_LO_genes, GO_list[GO_file], pval = 0.05)
    svglite(paste0("As4.1_Ren1-LO_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

Ren1_LO <- named_sct[, named_sct$named_cluster %in% "Ren1_LO"] 
table(Ren1_LO$orig.ident)
```

| cluster | ctrl (N) | H89 (N) |
|---------|----------|---------|
| Ren1_LO | 90       | 324     |

EMT is the highest enriched GO term in the Ren1_LO cluster (which is predominantly treated cells).

 - G2M/E2F/MYC targets all DOWN
 - EMT/P53/Myogenesis/TNFa/Apoptosis/Inflammatory response all UP

#### Cluster 4 == Sec24a|Arid3b+

Marked by the highest and near exclusive expression of *Sec24a* and *Arid3b*.

```R
cluster4_genes <- cluster_markers[cluster_markers$cluster=="Sec24a|Arid3b+",]$avg_log2FC
names(cluster4_genes) = cluster_markers[cluster_markers$cluster=="Sec24a|Arid3b+",]$gene
cluster4_genes = sort(cluster4_genes, decreasing = TRUE)
cluster4_genes = cluster4_genes[!duplicated(names(cluster4_genes))]
head(cluster4_genes)

for (GO_file in names(GO_list)) {
    res = GSEA(cluster4_genes, GO_list[GO_file], pval = 0.05)
    svglite(paste0("As4.1_Sec24a-Arid3b+_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

cluster4 <- named_sct[, named_sct$named_cluster %in% "Sec24a|Arid3b+"] 
table(cluster4$orig.ident)
```

| cluster        | ctrl (N) | H89 (N) |
|----------------|----------|---------|
| Sec24a|Arid3b+ | 193      | 178     |

#### Cluster 5 == Gbp7+

Marked by the highest and near exclusive expression of Gbp7.

```R
cluster5_genes <- cluster_markers[cluster_markers$cluster=="Gbp7+",]$avg_log2FC
names(cluster5_genes) = cluster_markers[cluster_markers$cluster=="Gbp7+",]$gene
cluster5_genes = sort(cluster5_genes, decreasing = TRUE)
cluster5_genes = cluster5_genes[!duplicated(names(cluster5_genes))]
head(cluster5_genes)

for (GO_file in names(GO_list)) {
    res = GSEA(cluster5_genes, GO_list[GO_file], pval = 0.05)
    svglite(paste0("As4.1_Gbp7+_", GO_file, ".svg"), 
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(res$Plot)
    dev.off()
}

cluster5 <- named_sct[, named_sct$named_cluster %in% "Gbp7+"] 
table(cluster5$orig.ident)
```

| cluster        | ctrl (N) | H89 (N) |
|----------------|----------|---------|
| Sec24a|Arid3b+ | 193      | 178     |

| cluster_condition     | N    |
|-----------------------|------|
| all+_ctrl             | 1511 |
| all+_h89              | 1755 |
| Gbp7+_ctrl            | 168  |
| Gbp7+_h89             | 148  |
| Khdrbs3\|Acaa1b+_ctrl | 588  |
| Khdrbs3\|Acaa1b+_h89  | 595  |
| remodel_ctrl          | 561  |
| remodel_h89           | 728  |
| Ren1_LO_ctrl          | 90   |
| Ren1_LO_h89           | 323  |
| Sec24a\|Arid3b+_ctrl  | 192  |
| Sec24a\|Arid3b+_h89   | 177  |