## 6. Perform velocity analysis

### 6a. Setup velocyto

```console
conda create velocyto numpy scipy cython numba matplotlib scikit-learn h5py click 
conda activate /project/gomezlab/data/.conda/envs/velocyto
pip install --user --upgrade --no-deps velocyto
pip install --user --upgrade --no-deps loompy
pip install --user --upgrade --no-deps numpy_groupies
```

### 6b. Run velocyto on control cells (i.e. DMSO-treated)

```console
velocyto run10x -m genomes/mouse/mm10_rmsk.gtf \
 As4-1-control_scRNAseq \
 genomes/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf
```

### 6c. Run velocyto on treated cells (i.e. H89-treated)

```console
velocyto run10x -m genomes/mouse/mm10_rmsk.gtf \
 As4-1-H89_scRNAseq \
 genomes/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf
```

### 6d. Load and organize velocyto data in R

```R
ctrl_loom <- "As4-1-control_scRNAseq/velocyto/As4-1-control_scRNAseq.loom"
trt_loom  <- "As4-1-H89_scRNAseq/velocyto/As4-1-H89_scRNAseq.loom"
```

Create and set up the control data.
1. Load velocyto file
2. Fix names
3. Keep matches

```R

ctrl_velo <- ReadVelocity(file = ctrl_loom)
ctrl_velo_seurat <- as.Seurat(ctrl_velo)

ctrl_velo_seurat <- RenameCells(
    ctrl_velo_seurat,
    new.names = gsub('As4-1-control_scRNAseq:', '', colnames(ctrl_velo_seurat))
)
ctrl_velo_seurat <- RenameCells(
    ctrl_velo_seurat,
    new.names = gsub('x', '', colnames(ctrl_velo_seurat))
)
Idents(ctrl_velo_seurat) <- gsub('As4-1-control_scRNAseq:', '',
                                 Idents(ctrl_velo_seurat))
Idents(ctrl_velo_seurat) <- gsub('x', '', Idents(ctrl_velo_seurat))

rm(fixed_names)
for (cellname in colnames(ctrl_velo_seurat)) {
    hit <- grep(cellname, colnames(named_sct))
    if (length(hit) == 1) {
        new_name <- colnames(named_sct)[hit]
    } else {
        new_name <- cellname
    }
    if (exists('fixed_names')) {
        fixed_names <- c(fixed_names, new_name)
    } else {
        fixed_names <- new_name
    }
}

ctrl_velo_seurat <- RenameCells(
    ctrl_velo_seurat,
    new.names = fixed_names
)
names_to_keep <- fixed_names[grep('-', fixed_names)]

ctrl_velo_seurat[["CellName"]] <- colnames(ctrl_velo_seurat)
ctrl_velo_seurat_keep <- subset(ctrl_velo_seurat,
                                subset = CellName %in% names_to_keep)

genes_in_use   <- rownames(GetAssayData(named_sct, assay = "RNA", slot = "data")) 
spliced_ctrl   <- GetAssayData(ctrl_velo_seurat_keep, assay = "spliced") %>%
    .[rownames(.) %in% genes_in_use,] %>% CreateAssayObject()
unspliced_ctrl <- GetAssayData(ctrl_velo_seurat_keep, assay = "unspliced") %>%
    .[rownames(.) %in% genes_in_use,] %>% CreateAssayObject()
ambiguous_ctrl <- GetAssayData(ctrl_velo_seurat_keep, assay = "ambiguous") %>%
    .[rownames(.) %in% genes_in_use,] %>% CreateAssayObject()
```

Repeat for the treated data.

```R
trt_velo  <- ReadVelocity(file = trt_loom)
trt_velo_seurat <- as.Seurat(trt_velo)

trt_velo_seurat <- RenameCells(
    trt_velo_seurat,
    new.names = gsub('As4-1-H89_scRNAseq:', '', colnames(trt_velo_seurat))
)
trt_velo_seurat <- RenameCells(
    trt_velo_seurat,
    new.names = gsub('x', '', colnames(trt_velo_seurat))
)
Idents(trt_velo_seurat) <- gsub('As4-1-H89_scRNAseq:', '',
                                Idents(trt_velo_seurat))
Idents(trt_velo_seurat) <- gsub('x', '', Idents(trt_velo_seurat))

rm(fixed_names)
for (cellname in colnames(trt_velo_seurat)) {
    hit <- grep(cellname, colnames(named_sct))
    if (length(hit) == 1) {
        new_name <- colnames(named_sct)[hit]
    } else {
        new_name <- cellname
    }
    if (exists('fixed_names')) {
        fixed_names <- c(fixed_names, new_name)
    } else {
        fixed_names <- new_name
    }
}

trt_velo_seurat <- RenameCells(
    trt_velo_seurat,
    new.names = fixed_names
)
names_to_keep <- fixed_names[grep('-', fixed_names)]

trt_velo_seurat[["CellName"]] <- colnames(trt_velo_seurat)
trt_velo_seurat_keep <- subset(trt_velo_seurat,
                               subset = CellName %in% names_to_keep)

genes_in_use  <- rownames(GetAssayData(named_sct, assay = "RNA", slot = "data")) 
spliced_trt   <- GetAssayData(trt_velo_seurat_keep, assay = "spliced") %>%
    .[rownames(.) %in% genes_in_use,] %>% CreateAssayObject()
unspliced_trt <- GetAssayData(trt_velo_seurat_keep, assay = "unspliced") %>%
    .[rownames(.) %in% genes_in_use,] %>% CreateAssayObject()
ambiguous_trt <- GetAssayData(trt_velo_seurat_keep, assay = "ambiguous") %>%
    .[rownames(.) %in% genes_in_use,] %>% CreateAssayObject()
```

Merge the Seurat objects, ensure names are correct, then grab the appropriate assays.

```R
velo_seurat <- merge(ctrl_velo_seurat_keep, y=trt_velo_seurat_keep,
                     add.cell.id=c('ctrl', 'trt'))
# Fix names
check_names <- colnames(velo_seurat)
check_names <- gsub('ctrl_', '', check_names)
check_names <- gsub('trt_', '', check_names)
duplicated_names <- check_names[duplicated(check_names)]
colnames(velo_seurat)[grep(duplicated_names, colnames(velo_seurat))]

paste0("trt_", duplicated_names)
names_to_keep <- names_to_keep[!(names_to_keep %in% duplicated_names)]
trt_velo_seurat_keep <- subset(trt_velo_seurat_keep,
                               subset = CellName %in% names_to_keep)

# Merge again
velo_seurat <- merge(ctrl_velo_seurat_keep, y=trt_velo_seurat_keep,
                     add.cell.id=c('ctrl', 'trt'))
velo_seurat <- RenameCells(
    velo_seurat,
    new.names = gsub('ctrl_', '', colnames(velo_seurat))
)
velo_seurat <- RenameCells(
    velo_seurat,
    new.names = gsub('trt_', '', colnames(velo_seurat))
)

# Now, extract the merged spliced/unspliced/ambiguous assays
spliced <- GetAssayData(velo_seurat, assay = "spliced") %>%
    .[rownames(.) %in% genes_in_use,] %>% CreateAssayObject()
unspliced <- GetAssayData(velo_seurat, assay = "unspliced") %>%
    .[rownames(.) %in% genes_in_use,] %>% CreateAssayObject()
ambiguous <- GetAssayData(velo_seurat, assay = "ambiguous") %>%
    .[rownames(.) %in% genes_in_use,] %>% CreateAssayObject()

# Subset the original named_sct object to only match cellnames from the velo_seurat object
named_sct[["CellName"]] <- colnames(named_sct)
named_sct_keep <- subset(named_sct, subset = CellName %in% colnames(velo_seurat))

# Then, add these to the named_sct object
named_sct_keep[["spliced"]]   <- spliced
named_sct_keep[["unspliced"]] <- unspliced
named_sct_keep[["ambiguous"]] <- ambiguous
```

### 6e. Run velocity analysis on merged data in R

```R
named_sct_keep <- RunVelocity(object = named_sct_keep, deltaT = 1,
                              kCells = 25, fit.quantile = 0.02)

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = named_sct_keep)))
names(x = ident.colors) <- levels(x = named_sct_keep)
cell.colors <- ident.colors[Idents(object = named_sct_keep)]
names(x = cell.colors) <- colnames(x = named_sct_keep)
svglite('velocity_plot.svg')
show.velocity.on.embedding.cor(emb = Embeddings(object = named_sct_keep,
    reduction = "umap"), vel = Tool(object = named_sct_keep, 
    slot = "RunVelocity"), n = 200, scale = "sqrt",
    cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE,
    min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1)
dev.off()
```

Plot marker genes on velocity embeddings

```R
emat <- named_sct_keep$spliced
nmat <- named_sct_keep$unspliced
cluster.label <- named_sct_keep$named_cluster
emat <- filter.genes.by.cluster.expression(emat, cluster.label,
    min.max.cluster.average = quantile(rowMeans(emat), 0.25))
nmat <- filter.genes.by.cluster.expression(nmat, cluster.label,
    min.max.cluster.average = quantile(rowMeans(nmat), 0.25))
fit.quantile <- 0.02
gene_signature_markers <- c("Gbp7", "Khdrbs3", "Acaa1b", "Sec24a", "Arid3b")
for (gene in gene_signature_markers) {
    svglite(paste0(gene, '_velocity_plot_kCells200.svg'), width=30,
            bg = "transparent", fix_text_size = FALSE, pointsize = 10)
    print(
        gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells = 200, 
            kGenes=1, fit.quantile=fit.quantile,
            cell.emb=Embeddings(object = named_sct_keep, reduction = "umap"),
            cell.colors=ac(x = cell.colors, alpha = 0.5),
            show.gene=gene,
            do.par=T)
         )
    dev.off()
}
```

Plot frequency distributions of cell populations.

```R
all_cells <- as.data.table(named_sct_keep$CellName)
all_cells[, cluster := named_sct_keep$named_cluster]
all_cells[, condition := named_sct_keep$orig.ident]

cell_distribution <- as.data.table(table(all_cells$condition,
                                         all_cells$cluster))
colnames(cell_distribution) <- c("condition", "cluster", "N")
cell_distribution$condition <- factor(cell_distribution$condition,
                                      levels = c("ctrl", "h89"))
cell_distribution_plot <- ggplot(cell_distribution,
    aes(x = cluster, y = N, fill = condition)) + 
  geom_bar(position = "dodge", stat = "identity", color="black") +
  scale_fill_manual(values = c("#000436", "#C732D5")) +
  labs(
        x = "Mouse kidney developmental timepoint",
        y = "Frequency"
      ) +
    custom_theme() +
    theme(legend.position = "right")
svglite(width=20,
        paste0("cell_distribution_frequency_side_barplot_",
               Sys.Date() ,".svg"),
        bg = "transparent", fix_text_size = FALSE, pointsize = 6)
print(cell_distribution_plot)
dev.off()
```