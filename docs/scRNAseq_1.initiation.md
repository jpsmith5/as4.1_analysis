## 1a. Initialize environment

Load necessary libraries and set initial paths to data files.

```R
library(dplyr)
library(patchwork)
library(svglite)
library(ggplot2)
library(sctransform)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(velocyto.R)
library(SeuratWrappers)
library(glmGamPoi)
library(data.table)
library(harmony)
library(fastSave)
library(ggstatsplot)
library(GSReg)
library(GSBenchMark)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ComplexHeatmap)
library(reshape2)
library(reticulate)
library(Rmagic)

set.seed(99)

sample_name   <- "As4.1_H89"
processed_dir <- "$HOME/processed/scRNA"

ctrl_matrix_path <- file.path(processed_dir,
    "As4-1-control_scRNAseq/outs/filtered_feature_bc_matrix/")
h89_matrix_path <- file.path(processed_dir,
    "As4-1-H89_scRNAseq/outs/filtered_feature_bc_matrix/")
```

## 1b. Set up a few helper functions

```R
custom_theme <- function(base_family = "sans", ...){
  theme_classic(base_family = base_family, base_size = 8, ...) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    aspect.ratio = 1,
    legend.position = "top",
    legend.justification="center",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )
}

get_random_distinct_colors <- function(ncolor, seed = 100) {
  require(uniformly)
  set.seed(seed)
  rgb_mat <- runif_in_cube(n=ncolor,d=3,O=rep(0.5,3),r=0.5)
  rgb(r=rgb_mat[,1],g=rgb_mat[,2],b=rgb_mat[,3])
}

get_random_grid_colors <- function(ncolor,seed = 100) {
  require(uniformly)
  set.seed(seed)
  ngrid <- ceiling(ncolor^(1/3))
  x <- seq(0,1,length=ngrid+1)[1:ngrid]
  dx <- (x[2] - x[1])/2
  x <- x + dx
  origins <- expand.grid(x,x,x)
  nbox <- nrow(origins) 
  RGB <- vector("numeric",nbox)
  for(i in seq_len(nbox)) {
    rgb <- runif_in_cube(n=1,d=3,O=as.numeric(origins[i,]),r=dx)
    RGB[i] <- rgb(rgb[1,1],rgb[1,2],rgb[1,3])
  }
  index <- sample(seq(1,nbox),ncolor)
  RGB[index]
}

pairwise_theme <- function(base_family = "sans", ...){
  theme_classic(base_family = base_family, base_size = 8, ...) +
  theme(
    axis.line = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    aspect.ratio = 1,
    legend.position = "top",
    legend.justification="center",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )
}
 
plotGEX <- function(seurat_object, gene, file_name, orig=FALSE, ...) {
    gene_plot  <- FeaturePlot(seurat_object,
                              features = paste0("sct_", gene),
                              reduction = 'umap',
                              cols = c("#D2D2D2", "#006CD1"),
                              raster=TRUE) +
                    theme_void() +
                    pairwise_theme() +
                    theme(axis.title.x = NULL)
    if (orig) {
        vln_plot   <- VlnPlot(seurat_object,
                              assay = "SCT",
                              cols = c("#A2A2A0", "#C30000"),
                              features = gene,
                              pt.size = 0.01,
                              ...) +
                        pairwise_theme() +
                        theme(aspect.ratio = NULL) +
                        NoLegend()
    } else {
        vln_plot   <- VlnPlot(seurat_object,
                              assay = "SCT",
                              cols = pal_named,
                              features = gene,
                              pt.size = 0.01,
                              ...) +
                        pairwise_theme() +
                        theme(aspect.ratio = NULL) +
                        NoLegend()
    }
    vln_plot2  <- VlnPlot(seurat_object,
                          assay = "SCT",
                          split.by = "orig.ident",
                          features = gene,
                          cols = c("#FFC20A", "#0C7BDC", "#A2A2A0"),
                          pt.size = 0.01,
                          ...) +
                    guides(color = guide_legend(direction = "horizontal")) +
                    pairwise_theme() +
                    theme(aspect.ratio = NULL)
    svglite(file_name, height=4, width=12, bg = "transparent",
            fix_text_size = FALSE, pointsize = 6)
    print(gene_plot | vln_plot | vln_plot2)
    dev.off()
}

GSEA <- function(gene_list, GO_file, pval) {
  set.seed(99)
  library(dplyr)
  library(fgsea)

  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list <- gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list <- sort(gene_list, decreasing = TRUE)
  }
  myGO <- fgsea::gmtPathways(GO_file)

  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=15,
                        maxSize=400,
                        nPermSimple = 10000,
                        nproc = 2) %>% 
                  as.data.frame() %>% 
                  dplyr::filter(padj < !!pval) %>% 
                  arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))

  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list)
  fgRes <- fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))

  fgRes$Enrichment <- ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes <- rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))

  total_up   <- sum(fgRes$Enrichment == "Up-regulated")
  total_down <- sum(fgRes$Enrichment == "Down-regulated")
  header     <- paste0("Top 10 (Total pathways: Up=", total_up,",
                        Down=", total_down, ")")

  colos <- setNames(c("firebrick2", "dodgerblue2"),
                    c("Up-regulated", "Down-regulated"))

  g1 <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
      geom_point( aes(fill = Enrichment, size = size), shape=21) +
      scale_fill_manual(values = colos ) +
      scale_size_continuous(range = c(2,10)) +
      geom_hline(yintercept = 0) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score", title=header) +
      custom_theme() +
      theme(legend.position = "right",
            panel.grid.major = element_line(color = "gray",
                                            linewidth = 0.1,
                                            linetype = 2))

  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}

removezeros <- function(p) {
    print(sum(apply(p,1,sd)==0,na.rm=T))
    keepIndex = (apply(p,1,sd)!=0)
    q = p[keepIndex,]
    print(dim(q))
    r <- t(q)
    return(r)
}

do_eva <- function(data, type, pathway, pheno) {
  ## transpose back for EVA
  edata <- t(data[[1]])
  exprsdata <- as.matrix(edata)
  rownames(exprsdata) <- rownames(edata)
  
  phenobin <- as.integer(as.factor(pheno))
  ## classic eva
  VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=exprsdata, pathways=pathway, phenotypes=factor(phenobin)) 
  pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);
  kendall <- lapply(VarAnKendallV, function (x) x[1:20])
  dd  <-  as.data.frame(matrix(unlist(kendall), nrow=length(unlist(kendall[1]))))
  rownames(dd) <- make.unique(names(kendall[[1]]))
  colnames(dd) <- names(kendall)
  dd <- as.data.frame(t(dd))
  dd$pathway_name <- rownames(dd)
  
  ## add meta information
  dd$metric <- "kendall_tau"
  dd$cell_line <- type
  dd$group1 <- "DMSO"
  dd$group2 <- "H89"
  return(dd)
}

plotEVAheatmap <- function(data, path) {
    ## split to join in long format
    E1 <- data[c(1,21,24)]
    E2 <- data[c(2,21,25)]
    colnames(E1) <- c("E","pathway","group")
    colnames(E2) <- c("E","pathway","group")
    plotdat <- rbind(E1,E2)

    mat  <- acast(plotdat, pathway ~ group , value.var='E')
    mats <- t(apply(mat,2,scale))
    rownames(mats) <- colnames(mat)
    colnames(mats) <- rownames(mat)

    ## add heatmap annotation
    ha_column = columnAnnotation(Treatment = as.factor(rownames(mats)))
    ha_row = rowAnnotation(Pathway = as.factor(colnames(mats)))
    m <- Heatmap(t(mats), col=viridis::inferno(5), name = "EVA statistics", 
            top_annotation=ha_column, 
            left_annotation=ha_row,
            border = TRUE,
            show_column_names = TRUE
    )

    return(m)
}

plotEVAboxplot <- function(data, path) {
    ## mult test correction
    data$adj.pval <- p.adjust(data$pvalue, method = "BH")
    print("number of significant pathways")
    print(nrow(data[data$adj.pval < 0.05,]))
    ## split to join in long format
    E1 <- data[c(1,21,24)]
    E2 <- data[c(2,21,25)]
    colnames(E1) <- c("E","pathway","group")
    colnames(E2) <- c("E","pathway","group")
    plotdat <- rbind(E1,E2)

    title <- paste("cell line", unique(data$cell_line), path)
    p <- ggplot(plotdat, aes(x=group, y=E, fill = group)) +
      geom_boxplot(outlier.shape = NA) + theme_bw() +
      theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
            axis.ticks.x=element_blank()) + ggtitle(title) +
      scale_fill_manual(values=c("blue", "red")) +
      geom_point(aes(fill = group, alpha = 0.8), size = 1.5, shape = 21, position = position_jitterdodge()) + 
      ylab("EVA statistic") + xlab("treatment")
    return(p)
}
```