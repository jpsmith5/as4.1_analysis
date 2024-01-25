## 3. Normalize and perform dimensionality reduction

```R
s.genes   <- stringr::str_to_sentence(cc.genes$s.genes)
g2m.genes <- stringr::str_to_sentence(cc.genes$g2m.genes)

ctrl_seurat <- CellCycleScoring(ctrl_seurat,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE)
h89_seurat <- CellCycleScoring(h89_seurat,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE)

ctrl <- SCTransform(ctrl_seurat, vst.flavor = "v2",
                    vars.to.regress = c("percent.mt", "percent.hemo",
                                        "S.Score", "G2M.Score"),
                    verbose = FALSE)
h89   <- SCTransform(h89_seurat, vst.flavor = "v2",
                     vars.to.regress = c("percent.mt", "percent.hemo",
                                         "S.Score", "G2M.Score"),
                     verbose = FALSE)

aoco_list <- list(ctrl = ctrl, treated = h89)
features  <- SelectIntegrationFeatures(object.list = aoco_list, nfeatures = 3000)
no_ribo   <- features[!grepl("^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa", features)]
no_riboORmito <- no_ribo[!grepl("^mt-", no_ribo)]
aoco_list <- PrepSCTIntegration(object.list = aoco_list,
                                anchor.features = no_riboORmito)
anchors <- FindIntegrationAnchors(object.list = aoco_list,
                                  normalization.method = "SCT",
                                  anchor.features = no_riboORmito)
combined_sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
combined_sct <- RunPCA(combined_sct, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.25, verbose = FALSE)

p1 <- DimPlot(combined_sct, reduction = "umap", group.by = "orig.ident") +
    custom_theme()
p2 <- DimPlot(combined_sct, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE) +
    custom_theme()
p3 <- DimPlot(combined_sct, reduction = "umap", group.by = "Phase",
              label = TRUE, repel = TRUE) +
    custom_theme()
p1 | p2 | p3

svglite(paste0(sample_name, "_UMAP_initial_", Sys.Date(), ".svg"),
        bg = "transparent", fix_text_size = FALSE, pointsize = 6)
p1 | p2 | p3
dev.off()

pal1 <- ArchR::paletteDiscrete(values = combined_sct$seurat_clusters)
pal2 <- get_random_grid_colors(length(unique(combined_sct$seurat_clusters)))
names(pal2) <- names(pal1)
pal <- pal2
```
