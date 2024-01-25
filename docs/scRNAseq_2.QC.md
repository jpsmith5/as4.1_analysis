## 2. Load and QC CellRanger processed data

### 2a. As4.1 Control (DMSO-treated) sample
```R
ctrl_RNA    <- Read10X(data.dir = file.path(ctrl_matrix_path))
ctrl_seurat <- CreateSeuratObject(counts = ctrl_RNA, project = "ctrl",
                                  min.cells = 3, min.features = 200)
ctrl_seurat[["percent.mt"]] <- PercentageFeatureSet(
    ctrl_seurat, pattern = "^mt-")
mtDNA_cutoff <- quantile(ctrl_seurat[["percent.mt"]]$percent.mt, 0.975)
    
ctrl_seurat[["percent.hemo"]] <- PercentageFeatureSet(object = ctrl_seurat,
                                                      pattern = "^Hb[^(p)]")
Hb_cutoff <- quantile(ctrl_seurat[["percent.hemo"]]$percent.hemo, 0.975)

plot1 <- FeatureScatter(ctrl_seurat,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.mt",
                        cols = 'gray')
svglite(paste0("ctrl", "_nCount-MT_", Sys.Date(), ".svg"))
plot1 + 
    geom_hline(yintercept = mtDNA_cutoff, col='red') + 
    custom_theme() +
    labs(x="RNA Count (N)", y="mtDNA (%)")
dev.off()

plot2 <- FeatureScatter(ctrl_seurat,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA",
                        cols = 'gray')
feature_upper_cutoff <- quantile(ctrl_seurat[["nFeature_RNA"]]$nFeature_RNA, 0.975)
feature_lower_cutoff <- quantile(ctrl_seurat[["nFeature_RNA"]]$nFeature_RNA, 0.025)
svglite(paste0("ctrl", "_nCount-nFeature_", Sys.Date(), ".svg"))
plot2 +
    geom_hline(yintercept = feature_upper_cutoff, col='red') + 
    geom_hline(yintercept = feature_lower_cutoff, col='red') + 
    custom_theme() +
    labs(x="RNA Count (N)", y="RNA Features (N)")
dev.off()

plot3 <- FeatureScatter(ctrl_seurat,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.hemo",
                        cols = 'gray')
svglite(paste0("ctrl", "_nCount-hemoglobin_", Sys.Date(), ".svg"))
plot3 + 
    geom_hline(yintercept = Hb_cutoff, col='red') + 
    custom_theme() +
    labs(x="RNA Count (N)", y="HbDNA (%)")
dev.off()

ctrl_seurat <- subset(ctrl_seurat,
                      subset = nFeature_RNA > feature_lower_cutoff &
                               nFeature_RNA < feature_upper_cutoff &
                               percent.mt   < mtDNA_cutoff &
                               percent.hemo < Hb_cutoff)
```

### 2b. As4.1 PKA inhibited (H89-treated) sample

```R
h89_RNA    <- Read10X(data.dir = file.path(h89_matrix_path))
h89_seurat <- CreateSeuratObject(counts = h89_RNA, project = "h89",
                                  min.cells = 3, min.features = 200)
h89_seurat[["percent.mt"]] <- PercentageFeatureSet(
    h89_seurat, pattern = "^mt-")
mtDNA_cutoff <- quantile(h89_seurat[["percent.mt"]]$percent.mt, 0.975)
    
h89_seurat[["percent.hemo"]] <- PercentageFeatureSet(object = h89_seurat,
                                                      pattern = "^Hb[^(p)]")
Hb_cutoff <- quantile(h89_seurat[["percent.hemo"]]$percent.hemo, 0.975)

plot1 <- FeatureScatter(h89_seurat,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.mt",
                        cols = 'gray')
svglite(paste0("h89", "_nCount-MT_", Sys.Date(), ".svg"))
plot1 + 
    geom_hline(yintercept = mtDNA_cutoff, col='red') + 
    custom_theme() +
    labs(x="RNA Count (N)", y="mtDNA (%)")
dev.off()

plot2 <- FeatureScatter(h89_seurat,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA",
                        cols = 'gray')
feature_upper_cutoff <- quantile(h89_seurat[["nFeature_RNA"]]$nFeature_RNA, 0.975)
feature_lower_cutoff <- quantile(h89_seurat[["nFeature_RNA"]]$nFeature_RNA, 0.025)
svglite(paste0("h89", "_nCount-nFeature_", Sys.Date(), ".svg"))
plot2 +
    geom_hline(yintercept = feature_upper_cutoff, col='red') + 
    geom_hline(yintercept = feature_lower_cutoff, col='red') + 
    custom_theme() +
    labs(x="RNA Count (N)", y="RNA Features (N)")
dev.off()

plot3 <- FeatureScatter(h89_seurat,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.hemo",
                        cols = 'gray')
svglite(paste0("h89", "_nCount-hemoglobin_", Sys.Date(), ".svg"))
plot3 + 
    geom_hline(yintercept = Hb_cutoff, col='red') + 
    custom_theme() +
    labs(x="RNA Count (N)", y="HbDNA (%)")
dev.off()

h89_seurat <- subset(h89_seurat,
                     subset = nFeature_RNA > feature_lower_cutoff &
                              nFeature_RNA < feature_upper_cutoff &
                              percent.mt   < mtDNA_cutoff &
                              percent.hemo < Hb_cutoff)
```
