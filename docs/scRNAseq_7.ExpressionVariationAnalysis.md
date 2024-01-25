## 7. Perform EVA analyses

### 7a. Extract matrices and impute

```R
ctrl_sct_keep <- named_sct_keep[, named_sct_keep$orig.ident == "ctrl"]
h89_sct_keep  <- named_sct_keep[, named_sct_keep$orig.ident == "h89"]
ctrl_GEX_keep <- GetAssayData(object = ctrl_sct_keep[["SCT"]], slot = "data")
h89_GEX_keep  <- GetAssayData(object = h89_sct_keep[["SCT"]], slot = "data")

ctrl_mat_keep <- as.matrix(ctrl_GEX_keep)
# add unique identifier to colnames
colnames(ctrl_mat_keep) <- paste(colnames(ctrl_mat_keep), "ctrl", sep = "_")

trt_mat_keep  <- as.matrix(h89_GEX_keep)
colnames(trt_mat_keep) <- paste(colnames(trt_mat_keep), "trt", sep = "_")

# pData for treatment
ctrl_pdat  <- data.frame("samples" = colnames(ctrl_GEX_keep),
                         "cell" = "As4.1", "treatment" = "DMSO")
trt_pdat   <- data.frame("samples" = colnames(h89_GEX_keep),
                         "cell" = "As4.1", "treatment" = "H89")
as4.1_pdat <- rbind(ctrl_pdat, trt_pdat)

exprsMat <- cbind(ctrl_mat_keep, trt_mat_keep)
save(exprsMat, file="exprsMat.Rda")

rownames(as4.1_pdat) <- as4.1_pdat$samples
save(as4.1_pdat, file="phenoData.Rda")

# generate fData
fdat <- toupper(rownames(ctrl_GEX_keep))
# set gene short names as rownames
rownames(fdat) <- toupper(rownames(ctrl_GEX_keep))
fdat <- as.data.frame(fdat)
fdat$feature <- fdat
fdat$feature <- "NA"
colnames(fdat) <- c("gene_short_name", "feature")
save(fdat, file="featureData.Rda")

as4.1_mat   <- removezeros(exprsMat)
m_as4.1_mat <- Rmagic::magic(as4.1_mat)
save(m_as4.1_mat, file = "magic_As4.1.rda")

# Label this matrix with cluster IDs in addition to CTRL vs H89
impute_mat <- m_as4.1_mat$result
impute_dt <- as.data.table(impute_mat, keep.rownames=TRUE)
tmp <- tstrsplit(impute_dt$rn, "_", fixed=TRUE)
impute_dt[, cell_name := paste0(tmp[[1]], "_", tmp[[2]])]
impute_dt[, treatment := tmp[[3]]]
```

### 7b. Save imputed matrix for each population of cells

```R
## all+
all_cellnames <- named_sct_keep[, named_sct_keep$named_cluster %in% "all+"]$CellName 
all_impute <- impute_dt[impute_dt$cell_name %in% all_cellnames, ]
all_impute[,cell_name:=NULL]
all_impute[,treatment:=NULL]
all_rownames <- all_impute$rn
all_impute[,rn:=NULL]

all_df <- as.data.frame(all_impute)
rownames(all_df) <- all_rownames
save(all_df, file = "magic_As4.1_all.rda")

## remodel+
remodel_cellnames <- named_sct_keep[, named_sct_keep$named_cluster %in% "remodel"]$CellName 
remodel_impute <- impute_dt[impute_dt$cell_name %in% remodel_cellnames, ]
remodel_impute[,cell_name:=NULL]
remodel_impute[,treatment:=NULL]
remodel_rownames <- remodel_impute$rn
remodel_impute[,rn:=NULL]

remodel_df <- as.data.frame(remodel_impute)
rownames(remodel_df) <- remodel_rownames
save(remodel_df, file = "magic_As4.1_remodel.rda")

## Ren1_LO
ren1_lo_cellnames <- named_sct_keep[, named_sct_keep$named_cluster %in% "Ren1_LO"]$CellName 
ren1_lo_impute <- impute_dt[impute_dt$cell_name %in% ren1_lo_cellnames, ]
ren1_lo_impute[,cell_name:=NULL]
ren1_lo_impute[,treatment:=NULL]
ren1_rownames <- ren1_lo_impute$rn
ren1_lo_impute[,rn:=NULL]

ren1_lo_df <- as.data.frame(ren1_lo_impute)
rownames(ren1_lo_df) <- ren1_rownames
save(ren1_lo_df, file = "magic_As4.1_Ren1_LO.rda")

## Sec24a|Arid3b
sec24a_arid3b_cellnames <- named_sct_keep[, named_sct_keep$named_cluster %in% "Sec24a|Arid3b+"]$CellName 
sec24a_arid3b_impute <- impute_dt[impute_dt$cell_name %in% sec24a_arid3b_cellnames, ]
sec24a_arid3b_impute[,cell_name:=NULL]
sec24a_arid3b_impute[,treatment:=NULL]
sec24a_arid3b_rownames <- sec24a_arid3b_impute$rn
sec24a_arid3b_impute[,rn:=NULL]

sec24a_arid3b_df <- as.data.frame(sec24a_arid3b_impute)
rownames(sec24a_arid3b_df) <- sec24a_arid3b_rownames
save(sec24a_arid3b_df, file = "magic_As4.1_Sec24a-Arid3b.rda")

## Khdrbs3|Acaa1b+
khdrbs3_acaa1b_cellnames <- named_sct_keep[, named_sct_keep$named_cluster %in% "Khdrbs3|Acaa1b+"]$CellName 
khdrbs3_acaa1b_impute <- impute_dt[impute_dt$cell_name %in% khdrbs3_acaa1b_cellnames, ]
khdrbs3_acaa1b_impute[,cell_name:=NULL]
khdrbs3_acaa1b_impute[,treatment:=NULL]
khdrbs3_acaa1b_rownames <- khdrbs3_acaa1b_impute$rn
khdrbs3_acaa1b_impute[,rn:=NULL]

khdrbs3_acaa1b_df <- as.data.frame(khdrbs3_acaa1b_impute)
rownames(khdrbs3_acaa1b_df) <- khdrbs3_acaa1b_rownames
save(khdrbs3_acaa1b_df, file = "magic_As4.1_Khdrbs3_Acaa1b.rda")

## Gbp7+
gbp7_cellnames <- named_sct_keep[, named_sct_keep$named_cluster %in% "Gbp7+"]$CellName 
gbp7_impute <- impute_dt[impute_dt$cell_name %in% gbp7_cellnames, ]
gbp7_impute[,cell_name:=NULL]
gbp7_impute[,treatment:=NULL]
gbp7_rownames <- gbp7_impute$rn
gbp7_impute[,rn:=NULL]

gbp7_df <- as.data.frame(gbp7_impute)
rownames(gbp7_df) <- gbp7_rownames
save(gbp7_df, file = "magic_As4.1_Gbp7.rda")
```

### 7c. For each cell population, run EVA analysis

Perform this process on the CTRL vs H89 samples within each cluster. EVA analysis takes ~17 hours to run. 

For example, using the Ren1_LO population:
```R
set.seed(99)

load("magic_As4.1_Ren1_LO.rda")

do_eva <- function(data, type, pathway, pheno) {
   ## transpose back for EVA
   edata <- t(data)
   exprsdata <- as.matrix(edata)
   rownames(exprsdata) <- rownames(edata)

   phenobin <- as.integer(as.factor(pheno))
   ## classic eva
   VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=exprsdata, pathways=pathway,
                                      phenotypes=factor(phenobin))
   pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);
   kendall <- lapply(VarAnKendallV, function (x) x[1:20])
   dd  <-  as.data.frame(matrix(unlist(kendall),
                         nrow=length(unlist(kendall[1]))))
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

## get binary vectors for each cell line
split1 <- strsplit(rownames(ren1_lo_df),split='_', fixed=TRUE)
c1     <- sapply(split1, "[", 3)

GO_hallmarks     <- "msigdb/mh.all.v2023.1.Mm.symbols.gmt"
GO_canonical     <- "msigdb/m2.all.v2023.1.Mm.symbols.gmt"
GO_TF_targets    <- "msigdb/m3.all.v2023.1.Mm.symbols.gmt"
GO_same_ontology <- "msigdb/m5.all.v2023.1.Mm.symbols.gmt"
GO_cell_sig      <- "msigdb/m8.all.v2023.1.Mm.symbols.gmt"

GO_list <- c(GO_hallmarks, GO_canonical, GO_TF_targets,
             GO_same_ontology, GO_cell_sig)
names(GO_list) <- c("hallmarks", "canonical", "TF_targets",
                    "same_ontology", "cell_signature")

library(GSReg)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
GO1 <- fgsea::gmtPathways(GO_list['hallmarks'])
Ren1_LO <- do_eva(data = ren1_lo_df, type = "As4.1", pathway = GO1, pheno = c1)
save(Ren1_LO, file = "EVA_Ren1_LO_hallmark.rda")

GO2 <- fgsea::gmtPathways(GO_list['canonical'])
Ren1_LO <- do_eva(data = ren1_lo_df, type = "As4.1", pathway = GO2, pheno = c1)
save(Ren1_LO, file = "EVA_Ren1_LO_canonical.rda")
```

### 7d. After evaluting ontologies for each population, combine and contrast.

```R
cluster_dfs <- c("gbp7_eva_hall", "gbp7_eva_canon", 
                 "khdrbs3_acaa1b_eva_hall", "khdrbs3_acaa1b_eva_canon",
                 "sec24a_arid3b_eva_hall", "sec24a_arid3b_eva_canon",
                 "remodel_eva_hall", "remodel_eva_canon",
                 "ren1_lo_eva_hall", "ren1_lo_eva_canon")
```

full dataset: H89 leads to slightly increased heterogeneity
 - all+:  H89 leads to increased heterogeneity
 - remodel: H89 leads to slightly increased heterogeneity
 - Ren1_LO: H89 has majorly reduced heterogeneity
 - Sec24a|Arid3b: H89 has reduced heterogeneity
 - Khdrbs3|Acaa1b: H89 leads to slightly increased heterogeneity
 - Gbp7+: H89 has reduced heterogeneity

Generate figures of EVA evaluations.

```R
for (file in cluster_dfs) {
    if (grepl('canon', file)) {
        pathname <- "canonical"
    } else {
        pathname <- "hallmark"
    }
    svglite(paste0("EVA_", file, "_boxplots.svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 6)
    print(plotEVAboxplot(data = get(file), path = pathname))
    dev.off()
    
    svglite(paste0("EVA_", file, "_heatmap.svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 6)
    print(plotEVAheatmap(data = get(file), path = pathname))
    dev.off()
}
```