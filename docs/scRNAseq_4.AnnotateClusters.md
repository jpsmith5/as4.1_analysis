## 4a. Identify cluster markers

```R
all_pos_markers <- FindAllMarkers(combined_sct, only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.25,
                                  random.seed=99)

top10_pos <- data.table(all_pos_markers %>% group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC))
top25_pos <- data.table(all_pos_markers %>% group_by(cluster) %>%
    top_n(n = 25, wt = avg_log2FC))
top50_pos <- data.table(all_pos_markers %>% group_by(cluster) %>%
    top_n(n = 50, wt = avg_log2FC))

fwrite(top10_pos,
    paste0(sample_name, "_Top10PositiveMarkers_", Sys.Date(), ".csv"),
    sep=",")
fwrite(top25_pos,
    paste0(sample_name, "_Top25PositiveMarkers_", Sys.Date(), ".csv"),
    sep=",")
fwrite(top50_pos,
    paste0(sample_name, "_Top50PositiveMarkers_", Sys.Date(), ".csv"),
    sep=",")

Idents(combined_sct) <- combined_sct@meta.data$'orig.ident'
combined_sct <- PrepSCTFindMarkers(combined_sct, assay="SCT")
orig_markers <- FindAllMarkers(combined_sct, only.pos = FALSE, assay="SCT",
                               min.pct = 0.25, logfc.threshold = 0.25,
                               test.use="MAST", random.seed=99)

ctrl_gene_list = orig_markers[orig_markers$cluster=="ctrl",]$avg_log2FC
names(ctrl_gene_list) = orig_markers[orig_markers$cluster=="ctrl",]$gene
ctrl_gene_list = sort(ctrl_gene_list, decreasing = TRUE)
ctrl_gene_list = ctrl_gene_list[!duplicated(names(ctrl_gene_list))]

h89_gene_list = orig_markers[orig_markers$cluster=="h89",]$avg_log2FC
names(h89_gene_list) = orig_markers[orig_markers$cluster=="h89",]$gene
h89_gene_list = sort(h89_gene_list, decreasing = TRUE)
h89_gene_list = h89_gene_list[!duplicated(names(h89_gene_list))]

fwrite(orig_markers[orig_markers$cluster=="ctrl",], 
       paste0(sample_name, "_CTRL_markers.csv"), sep=",")
fwrite(orig_markers[orig_markers$cluster=="h89",], 
       paste0(sample_name, "_TRT_markers.csv"), sep=",")
```

Evaluate renin lineage and smooth muscle cell markers

```R
renin_lineage_markers <- c("Foxd1", "Ren1", "Akr1b7")
for (gene in renin_lineage_markers) {
    svglite(paste0(sample_name, "_", gene, "_VlnPlot_", Sys.Date(), ".svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 6)
    print(VlnPlot(combined_sct, assay="SCT",
            features = gene, group.by="seurat_clusters", cols=pal) +
        custom_theme() +
        theme(legend.position="none"))
    dev.off()
}
smc_markers <- c("Acta2", "Myh11", "Tagln", "Cnn1", "Smtn")
for (gene in smc_markers) {
    svglite(paste0(sample_name, "_", gene, "_VlnPlot_", Sys.Date(), ".svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 6)
    print(VlnPlot(combined_sct, assay="SCT",
            features = gene, group.by="seurat_clusters", cols=pal) +
        custom_theme() +
        theme(legend.position="none"))
    dev.off()
}
```

 - Cluster 3 is predominantly treated samples, with correspondingly lower *Ren1* expression/

```R
svglite(paste0(sample_name, "_TopMarkers_VlnPlot_", Sys.Date(), ".svg"),
        bg = "transparent", fix_text_size = FALSE, pointsize = 4,
        height=24, width=24)
VlnPlot(combined_sct, features = c("Eprs",
                                   "Hist1h2ap",
                                   "Acaa1b",
                                   "Cd109",
                                   "Pttg1",
                                   "Gbp7"),
        assay="SCT",
        group.by="seurat_clusters",
        pt.size=1/32)
dev.off()

rank2markers <- c("Eif4ebp1", "Hist1h2ap", "Khdrbs3", "Dapk1",
                  "Cyb5r3", "Ccl4")

svglite(paste0(sample_name, "_rank2markers_VlnPlot_", Sys.Date(), ".svg"),
        bg = "transparent", fix_text_size = FALSE, pointsize = 4,
        height=24, width=24)
VlnPlot(combined_sct,
        features = rank2markers,
        assay="SCT",
        group.by="seurat_clusters",
        pt.size=1/32)
dev.off()

cluster0 <- c("Eif4ebp1", "Hspa9", "Eprs", "Papss1", "Cyb5r1", "Eif2s2",
              "1110038B12Rik", "Phgdh", "Ptgfr", "Tars")
VlnPlot(combined_sct,
        features = cluster0,
        assay="SCT",
        group.by="seurat_clusters",
        pt.size=1/32)

cluster4 <- c("Cyb5r3", "Sec24a", "Dennd2a", "Hmmr", "Arid3b")
cluster4 <- c("Prkch", "Emp1", "Pttg1", "Pkdcc", "Amotl2")

VlnPlot(combined_sct,
        features = cluster4,
        assay="SCT",
        group.by="seurat_clusters",
        pt.size=1/32)

```

 - Cluster 5: Top marker [Gbp7](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GBP7)
 
```R

```

## 4b. Annotate clusters

| cluster | abbreviation    | notes                                       |
|---------|-----------------|---------------------------------------------|
| 0       | all+            | no distinguishing gene enrichment           |
| 1       | remodel         | chromatin remodeling pathways               |
| 2       | Khdrbs3|Acaa1b+ | highest expression for *Khdrbs3* & *Acaa1b* |
| 3       | Ren1_LO         | mostly treated cells with low *Ren1*        |
| 4       | Sec24a|Arid3b+  | highest expression for *Sec24a* & *Arid3b*  |
| 5       | Gbp7+           | only cluster enriched for *Gbp7*            |


| cluster_condition | N    | marker_genes                  |
|-------------------|------|-------------------------------|
| all+              | 3266 | --                            |
| remodel           | 1289 | Hist1h2ap, Hist1h1b, Hist1h1e |
| Khdrbs3\|Acaa1b+  | 1183 | Khdrbs3, Acaa1b               |
| Ren1_LO           | 413  | Ren1                          |
| Sec24a\|Arid3b+   | 369  | Sec24a, Arid3b                |
| Gbp7+             | 316  | Gbp7                          |

```R
Idents(combined_sct) <- combined_sct$seurat_clusters
new_cluster_ids <- c("all+", "remodel", "Khdrbs3|Acaa1b+",
                     "Ren1_LO", "Sec24a|Arid3b+", "Gbp7+")
names(new_cluster_ids) <- c("0", "1", "2", "3", "4", "5")
named_sct <- copy(combined_sct)
named_sct <- RenameIdents(named_sct, new_cluster_ids)
named_sct[["named_cluster"]] <- Idents(named_sct)
table(named_sct$named_cluster)
```

| cluster | abbreviation    | N    |
|---------|-----------------|------|
| 0       | all+            | 3268 |
| 1       | remodel         | 1290 |
| 2       | Khdrbs3|Acaa1b+ | 1183 |
| 3       | Ren1_LO         | 414  |
| 4       | Sec24a|Arid3b+  | 371  |
| 5       | Gbp7+           | 318  |

Extract gene expression matrices

```R
ctrl_sct <- combined_sct[, combined_sct$orig.ident == "ctrl"]
h89_sct  <- combined_sct[, combined_sct$orig.ident == "h89"]
ctrl_GEX <- GetAssayData(object = ctrl_sct[["SCT"]], slot = "data")
h89_GEX  <- GetAssayData(object = h89_sct[["SCT"]], slot = "data")

ctrl_mat <- as.data.table(t(as.matrix(ctrl_GEX)), keep.rownames=TRUE)
ctrl_mat[, group:="ctrl"]

trt_mat <- as.data.table(t(as.matrix(h89_GEX)), keep.rownames=TRUE)
trt_mat[, group:="trt"]

# Grab index values and names for each cell annotation
rm(cell_indicies)
for (celltype in unique(named_sct$named_cluster)) {
    celltype_name <- celltype
    if (grepl("+", celltype, fixed=T)) {
        celltype = sub("+", "\\+", celltype, fixed=T)
    }
    indices   <- grep(paste0("^", celltype, "$"), Idents(named_sct))
    cellnames <- names(Idents(named_sct)[indices])
    tmp_dt <- data.table(cellname = cellnames,
                         index = indices,
                         named_cluster = rep(celltype_name, length(cellnames)))
    if (exists("cell_indicies")) {
        cell_indicies <- rbind(cell_indicies, tmp_dt)
    } else {
        cell_indicies <- tmp_dt
    }
}

ctrl_mat <- merge(ctrl_mat, cell_indicies, by.x="rn", by.y="cellname")
trt_mat  <- merge(trt_mat, cell_indicies, by.x="rn", by.y="cellname")

GEX_dt <- rbind(ctrl_mat, trt_mat)

fwrite(GEX_dt, "As4.1_H89-DMSO_GEX.csv", sep=",")
```

Plot marker genes

```R
cluster_marker_gene_list <- c("Ren1", "Hist1h2ap", "Hist1h1b", "Hist1h1e",
                              "Khdrbs3", "Acaa1b", "Sec24a", "Arid3b", "Gbp7")

g1 <- ggstatsplot::ggbetweenstats(
      data  = GEX_dt,
      x     = named_cluster,
      y     = Ren1,
      title = "Ren1",
      ggtheme = custom_theme() + theme(panel.grid.major = element_line(
                                        color = "gray",
                                        linewidth = 0.1,
                                        linetype = 2))
    )
svglite(paste0("Ren1_SCT_ggbetweenstats_", Sys.Date(), ".svg"),
        bg = "transparent", fix_text_size = FALSE, pointsize = 4)
print(g1)
dev.off()

gene_signature_colors <- c("#95CF88", "#C90076", "#873C3E",
                           "#2ED8C4", "#272146", "#2433E1")
names(gene_signature_colors) <- c("all+", "remodel", "Ren1_LO",
                                  "Khdrbs3|Acaa1b+", "Sec24a|Arid3b+", "Gbp7+")

for (gene in cluster_marker_gene_list) {
    base_plot <- ggplot(GEX_dt, aes(x=named_cluster, y=!!rlang::sym(gene)))
    svglite(paste0(gene, "_SCT_boxplot_", Sys.Date(), ".svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(base_plot + 
        geom_jitter(width = 0.125, stroke = 0, size=3,
                    aes(col=GEX_dt$named_cluster, group=GEX_dt$named_cluster, alpha=0.1)) +
        geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
        #stat_boxplot(geom ='errorbar', width = 0.5) +
        geom_violin(width = 0.50, fill = "transparent",
                    draw_quantiles = 0.5, col='#333333') +
        stat_summary(fun = "mean", geom = "point",
                     shape = 1, size = 2) +
        labs(x="", y="normalized counts", title=gene) +
        scale_color_manual(values=gene_signature_colors) +
        scale_fill_manual(values=gene_signature_colors) +
        custom_theme() +
        ylim(0,8))
    dev.off()
}
```