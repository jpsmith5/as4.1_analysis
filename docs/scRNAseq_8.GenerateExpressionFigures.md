## 8. Generate figures

From insight gained from companion ATAC-seq data, we can plot expression values for genes of interest in As4.1 cells.

### 8a. Generate expression plots for transcription factor (TF) family members

For AP-1 TFs:

```R
AP1_members <- c("Junb", "Nfe2l1", "Fos", "Jund", "Jdp2", "Nfe2l2", "Fosl1",
                 "Fosl2", "Mafg", "Maff", "Atf3", "Bach1", "Maf", "Mafk",
                 "Fosb", "Mafb", "Batf3", "Batf2", "Batf", "Mafa", "Bach2",
                 "Nfe2l3", "Nfe2", "Nrl")

# Individual plots for each gene across the cell populations.
for (gene in AP1_members) {
    base_plot <- ggplot(GEX_dt, aes(x=named_cluster, y=!!rlang::sym(gene)))
    svglite(paste0(gene, "_SCT_boxplot_", Sys.Date(), ".svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(base_plot + 
        geom_jitter(width = 0.125, stroke = 0, size=3,
                    aes(col=GEX_dt$named_cluster, group=GEX_dt$named_cluster, alpha=0.1)) +
        geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
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

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
AP1_present <- character()
for (gene in AP1_members) {
    if (length(GEX_dt[[gene]]) != 0) {
        AP1_present <- c(AP1_present, gene)
    }
}

GEX_AP1 <- subset(GEX_dt, select=c(AP1_present, "group", "named_cluster"))
GEX_AP1_reshape <- reshape2::melt(GEX_AP1)

AP1_summary <- as.data.table(GEX_AP1_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value),
          Median=median(value), Std=sd(value)))
levels <- AP1_summary[order(Mean, decreasing=T),]$variable
GEX_AP1_reshape$variable <- factor(GEX_AP1_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_AP1_reshape[GEX_AP1_reshape$group=="ctrl",],
                    aes(x=variable, y=value)) 
svglite(paste0("AP1_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="AP1 TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
# Plot the genes all together on single plot just within the "Ren1_LO" cluster
Ren1_LO_GEX <- GEX_dt[GEX_dt$named_cluster == "Ren1_LO", ]

AP1_present <- character()
for (gene in AP1_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        AP1_present <- c(AP1_present, gene)
    }
}

Ren1_LO_AP1 <- subset(Ren1_LO_GEX, select=c(AP1_present, "group", "named_cluster"))
Ren1_LO_AP1_reshape <- reshape2::melt(Ren1_LO_AP1)

AP1_summary <- as.data.table(Ren1_LO_AP1_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- AP1_summary[order(Mean, decreasing=T),]$variable
Ren1_LO_AP1_reshape$variable <- factor(Ren1_LO_AP1_reshape$variable, levels=levels)

base_plot <- ggplot(Ren1_LO_AP1_reshape, aes(x=variable, y=value)) 
svglite(paste0("AP1_SCT_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3, aes(col=gene_signature_colors[["Ren1_LO"]], alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="AP1 TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

For KLF TFs:

```R
KLF_members <- c("Klf1", "Klf2", "Klf3", "Klf4", "Klf5", "Klf6", "Klf7",
                 "Klf8", "Klf9", "Klf10", "Klf11", "Klf12", "Klf13", "Klf14",
                 "Klf15", "Klf16", "Klf17")
```

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
KLF_present <- character()
for (gene in KLF_members) {
    if (length(GEX_dt[[gene]]) != 0) {
        KLF_present <- c(KLF_present, gene)
    }
}

GEX_KLF <- subset(GEX_dt, select=c(KLF_present, "group", "named_cluster"))
GEX_KLF_reshape <- reshape2::melt(GEX_KLF)

KLF_summary <- as.data.table(GEX_KLF_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- KLF_summary[order(Mean, decreasing=T),]$variable
GEX_KLF_reshape$variable <- factor(GEX_KLF_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_KLF_reshape[GEX_KLF_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("KLF_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="KLF TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
KLF_present <- character()
for (gene in KLF_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        KLF_present <- c(KLF_present, gene)
    }
}

Ren1_LO_KLF <- subset(Ren1_LO_GEX, select=c(KLF_present, "group", "named_cluster"))
Ren1_LO_KLF_reshape <- reshape2::melt(Ren1_LO_KLF)

KLF_summary <- as.data.table(Ren1_LO_KLF_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- KLF_summary[order(Mean, decreasing=T),]$variable
Ren1_LO_KLF_reshape$variable <- factor(Ren1_LO_KLF_reshape$variable, levels=levels)

base_plot <- ggplot(Ren1_LO_KLF_reshape, aes(x=variable, y=value)) 
svglite(paste0("KLF_SCT_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3, aes(col=gene_signature_colors[["Ren1_LO"]], alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="KLF TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

For SP TFs:

```R
SP_members <- c("Sp1", "Sp2", "Sp3", "Sp4", "Sp5", "Sp6", "Sp7", "Sp8", "Sp9")
```

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
SP_present <- character()
for (gene in SP_members) {
    if (length(GEX_dt[[gene]]) != 0) {
        SP_present <- c(SP_present, gene)
    }
}

GEX_SP <- subset(GEX_dt, select=c(SP_present, "group", "named_cluster"))
GEX_SP_reshape <- reshape2::melt(GEX_SP)

SP_summary <- as.data.table(GEX_SP_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- SP_summary[order(Mean, decreasing=T),]$variable
GEX_SP_reshape$variable <- factor(GEX_SP_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_SP_reshape[GEX_SP_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("SP_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="SP TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
SP_present <- character()
for (gene in SP_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        SP_present <- c(SP_present, gene)
    }
}

Ren1_LO_SP <- subset(Ren1_LO_GEX, select=c(SP_present, "group", "named_cluster"))
Ren1_LO_SP_reshape <- reshape2::melt(Ren1_LO_SP)

SP_summary <- as.data.table(Ren1_LO_SP_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- SP_summary[order(Mean, decreasing=T),]$variable
Ren1_LO_SP_reshape$variable <- factor(Ren1_LO_SP_reshape$variable, levels=levels)

base_plot <- ggplot(Ren1_LO_SP_reshape, aes(x=variable, y=value)) 
svglite(paste0("SP_SCT_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3, aes(col=gene_signature_colors[["Ren1_LO"]], alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="SP TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

For ZBTB TFs:

```R
ZBTB_members <- c("Patz1", "Hic1", "Bcl6", "Bcl6b", "Plzf", "Znf131", 
                  "Zbtb1", "Zbtb2", "Zbtb3", "Zbtb4", "Zbtb5", "Zbtb6",
                  "Zbtb7a", "Zbtb7b", "Zbtb11", "Zbtb12", "Zbtb13",
                  "Zbtb14", "Zbtb15", "Zbtb17", "Zbtb18", 
                  "Zbtb20", "Zbtb21", "Zbtb24", "Zbtb25", "Zbtb26", "Zbtb27",
                  "Zbtb28", "Zbtb29", "Zbtb30", "Zbtb31", "Zbtb32", "Zbtb33",
                  "Zbtb34", "Zbtb36", "Zbtb37", "Zbtb38", "Zbtb39", "Zbtb40")
```

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
ZBTB_present <- character()
for (gene in ZBTB_members) {
    if (length(GEX_dt[[gene]]) != 0) {
        ZBTB_present <- c(ZBTB_present, gene)
    }
}

GEX_ZBTB <- subset(GEX_dt, select=c(ZBTB_present, "group", "named_cluster"))
GEX_ZBTB_reshape <- reshape2::melt(GEX_ZBTB)
# Sort by median expression
ZBTB_summary <- as.data.table(GEX_ZBTB_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))

levels <- ZBTB_summary[order(Mean, decreasing=T),]$variable

GEX_ZBTB_reshape$variable <- factor(GEX_ZBTB_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_ZBTB_reshape[GEX_ZBTB_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("ZBTB_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="ZBTB TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
ZBTB_present <- character()
for (gene in ZBTB_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        ZBTB_present <- c(ZBTB_present, gene)
    }
}

Ren1_LO_ZBTB <- subset(Ren1_LO_GEX, select=c(ZBTB_present, "group", "named_cluster"))
Ren1_LO_ZBTB_reshape <- reshape2::melt(Ren1_LO_ZBTB)
# Sort by median expression
ZBTB_summary <- as.data.table(Ren1_LO_ZBTB_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))

levels <- ZBTB_summary[order(Mean, decreasing=T),]$variable

Ren1_LO_ZBTB_reshape$variable <- factor(Ren1_LO_ZBTB_reshape$variable, levels=levels)

base_plot <- ggplot(Ren1_LO_ZBTB_reshape, aes(x=variable, y=value)) 
svglite(paste0("ZBTB_SCT_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3, aes(col=gene_signature_colors[["Ren1_LO"]], alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="ZBTB TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

For ETS TFs:

```R
ETS_members <- c("Elf1", "Elf2", "Elf4", "Gabpa", "Erg", "Fli1", "Fev", "Erf",
                 "Etv3", "Elf3", "Elf5", "Ese3", "Ets1", "Ets2", "Spdef",
                 "Etv4", "Etv5", "Etv1", "Etv2", "Spi1", "Spib", "Elk3", "Etv6",
                 "Etv7")
```

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
ETS_present <- character()
for (gene in ETS_members) {
    if (length(GEX_dt[[gene]]) != 0) {
        ETS_present <- c(ETS_present, gene)
    }
}

GEX_ETS <- subset(GEX_dt, select=c(ETS_present, "group", "named_cluster"))
GEX_ETS_reshape <- reshape2::melt(GEX_ETS)
# Sort by mean expression
ETS_summary <- as.data.table(GEX_ETS_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))

levels <- ETS_summary[order(Mean, decreasing=T),]$variable

GEX_ETS_reshape$variable <- factor(GEX_ETS_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_ETS_reshape[GEX_ETS_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("ETS_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="ETS TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
ETS_present <- character()
for (gene in ETS_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        ETS_present <- c(ETS_present, gene)
    }
}

Ren1_LO_ETS <- subset(Ren1_LO_GEX, select=c(ETS_present, "group", "named_cluster"))
Ren1_LO_ETS_reshape <- reshape2::melt(Ren1_LO_ETS)
ETS_summary <- as.data.table(Ren1_LO_ETS_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))

levels <- ETS_summary[order(Mean, decreasing=T),]$variable

Ren1_LO_ETS_reshape$variable <- factor(Ren1_LO_ETS_reshape$variable, levels=levels)

base_plot <- ggplot(Ren1_LO_ETS_reshape, aes(x=variable, y=value)) 
svglite(paste0("ETS_SCT_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3, aes(col=gene_signature_colors[["Ren1_LO"]], alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="ETS TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

For Maz TFs:

```R
MAZ_members <- c("Maz")
MAZ_present <- character()
for (gene in MAZ_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        MAZ_present <- c(MAZ_present, gene)
    }
}
```

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
GEX_MAZ <- subset(GEX_dt, select=c(MAZ_present, "group", "named_cluster"))
GEX_MAZ_reshape <- reshape2::melt(GEX_MAZ)
MAZ_summary <- as.data.table(GEX_MAZ_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))

levels <- MAZ_summary[order(Mean, decreasing=T),]$variable

GEX_MAZ_reshape$variable <- factor(GEX_MAZ_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_MAZ_reshape[GEX_MAZ_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("MAZ_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="MAZ TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
MAZ_present <- character()
for (gene in MAZ_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        MAZ_present <- c(MAZ_present, gene)
    }
}

Ren1_LO_MAZ <- subset(Ren1_LO_GEX, select=c(MAZ_present, "group", "named_cluster"))
Ren1_LO_MAZ_reshape <- reshape2::melt(Ren1_LO_MAZ)
MAZ_summary <- as.data.table(Ren1_LO_MAZ_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))

levels <- MAZ_summary[order(Mean, decreasing=T),]$variable

Ren1_LO_MAZ_reshape$variable <- factor(Ren1_LO_MAZ_reshape$variable, levels=levels)


base_plot <- ggplot(Ren1_LO_MAZ_reshape[Ren1_LO_MAZ_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("MAZ_SCT_CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="MAZ TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

For P53 TFs:

```R
P53_members <- c("Trp63", "Trp53", "Trp73")
```

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
P53_present <- character()
for (gene in P53_members) {
    if (length(GEX_dt[[gene]]) != 0) {
        P53_present <- c(P53_present, gene)
    }
}

GEX_P53 <- subset(GEX_dt, select=c(P53_present, "group", "named_cluster"))
GEX_P53_reshape <- reshape2::melt(GEX_P53)

P53_summary <- as.data.table(GEX_P53_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- P53_summary[order(Mean, decreasing=T),]$variable
GEX_P53_reshape$variable <- factor(GEX_P53_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_P53_reshape[GEX_P53_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("P53_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="P53 TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
P53_present <- character()
for (gene in P53_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        P53_present <- c(P53_present, gene)
    }
}

Ren1_LO_P53 <- subset(Ren1_LO_GEX, select=c(P53_present, "group", "named_cluster"))
Ren1_LO_P53_reshape <- reshape2::melt(Ren1_LO_P53)

P53_summary <- as.data.table(Ren1_LO_P53_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- P53_summary[order(Mean, decreasing=T),]$variable
Ren1_LO_P53_reshape$variable <- factor(Ren1_LO_P53_reshape$variable, levels=levels)

base_plot <- ggplot(Ren1_LO_P53_reshape[Ren1_LO_P53_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("P53_SCT_CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="P53 TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

For TEAD TFs:

```R
TEAD_members <- c("Tead4", "Tead1", "Tead2", "Tead3")
```

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
TEAD_present <- character()
for (gene in TEAD_members) {
    if (length(GEX_dt[[gene]]) != 0) {
        TEAD_present <- c(TEAD_present, gene)
    }
}

GEX_TEAD <- subset(GEX_dt, select=c(TEAD_present, "group", "named_cluster"))
GEX_TEAD_reshape <- reshape2::melt(GEX_TEAD)

TEAD_summary <- as.data.table(GEX_TEAD_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- TEAD_summary[order(Mean, decreasing=T),]$variable
GEX_TEAD_reshape$variable <- factor(GEX_TEAD_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_TEAD_reshape[GEX_TEAD_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("TEAD_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="TEAD TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
TEAD_present <- character()
for (gene in TEAD_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        TEAD_present <- c(TEAD_present, gene)
    }
}

Ren1_LO_TEAD <- subset(Ren1_LO_GEX, select=c(TEAD_present, "group", "named_cluster"))
Ren1_LO_TEAD_reshape <- reshape2::melt(Ren1_LO_TEAD)

TEAD_summary <- as.data.table(Ren1_LO_TEAD_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- TEAD_summary[order(Mean, decreasing=T),]$variable
Ren1_LO_TEAD_reshape$variable <- factor(Ren1_LO_TEAD_reshape$variable, levels=levels)

base_plot <- ggplot(Ren1_LO_TEAD_reshape[Ren1_LO_TEAD_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("TEAD_SCT_CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="TEAD TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

For RUNX TFs:

```R
RUNX_members <- c("Runx1", "Runx2", "Runx3")
```

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
RUNX_present <- character()
for (gene in RUNX_members) {
    if (length(GEX_dt[[gene]]) != 0) {
        RUNX_present <- c(RUNX_present, gene)
    }
}

GEX_RUNX <- subset(GEX_dt, select=c(RUNX_present, "group", "named_cluster"))
GEX_RUNX_reshape <- reshape2::melt(GEX_RUNX)

RUNX_summary <- as.data.table(GEX_RUNX_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- RUNX_summary[order(Mean, decreasing=T),]$variable
GEX_RUNX_reshape$variable <- factor(GEX_RUNX_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_RUNX_reshape[GEX_RUNX_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("RUNX_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="RUNX TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
RUNX_present <- character()
for (gene in RUNX_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        RUNX_present <- c(RUNX_present, gene)
    }
}

Ren1_LO_RUNX <- subset(Ren1_LO_GEX, select=c(RUNX_present, "group", "named_cluster"))
Ren1_LO_RUNX_reshape <- reshape2::melt(Ren1_LO_RUNX)

RUNX_summary <- as.data.table(Ren1_LO_RUNX_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- RUNX_summary[order(Mean, decreasing=T),]$variable
Ren1_LO_RUNX_reshape$variable <- factor(Ren1_LO_RUNX_reshape$variable, levels=levels)

base_plot <- ggplot(Ren1_LO_RUNX_reshape[Ren1_LO_RUNX_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("RUNX_SCT_CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="RUNX TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

For ROR TFs:

```R
ROR_members <- c("Rora", "Rorb", "Rorc")
```

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
ROR_present <- character()
for (gene in ROR_members) {
    if (length(GEX_dt[[gene]]) != 0) {
        ROR_present <- c(ROR_present, gene)
    }
}

GEX_ROR <- subset(GEX_dt, select=c(ROR_present, "group", "named_cluster"))
GEX_ROR_reshape <- reshape2::melt(GEX_ROR)

ROR_summary <- as.data.table(GEX_ROR_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- ROR_summary[order(Mean, decreasing=T),]$variable
GEX_ROR_reshape$variable <- factor(GEX_ROR_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_ROR_reshape[GEX_ROR_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("ROR_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="ROR TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
ROR_present <- character()
for (gene in ROR_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        ROR_present <- c(ROR_present, gene)
    }
}

Ren1_LO_ROR <- subset(Ren1_LO_GEX, select=c(ROR_present, "group", "named_cluster"))
Ren1_LO_ROR_reshape <- reshape2::melt(Ren1_LO_ROR)

ROR_summary <- as.data.table(Ren1_LO_ROR_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- ROR_summary[order(Mean, decreasing=T),]$variable
Ren1_LO_ROR_reshape$variable <- factor(Ren1_LO_ROR_reshape$variable, levels=levels)

base_plot <- ggplot(Ren1_LO_ROR_reshape[Ren1_LO_ROR_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("ROR_SCT_CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="ROR TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

For ZNF TFs:

```R
ZNF_members <- fread('ZNF_gene_list.csv', header=F))
ZNF_members <- ZNF_members$V1

ZNF_members <- colnames(GEX_dt)[grep('Zfp', colnames(GEX_dt))]

# In all scRNA-seq from As4.1
ZNF_present <- character()
for (gene in ZNF_members) {
    if (length(GEX_dt[[gene]]) != 0) {
        ZNF_present <- c(ZNF_present, gene)
    }
}

```

In all scRNA-seq from control (DMSO-treated) As4.1 cells.

```R
GEX_ZNF <- subset(GEX_dt, select=c(ZNF_present, "group", "named_cluster"))
GEX_ZNF_reshape <- reshape2::melt(GEX_ZNF)

ZNF_summary <- as.data.table(GEX_ZNF_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))

exp_list <- c(as.character(ZNF_summary[Mean>0 & Median>1,]$variable), "Zfp422", "Zfp467", "Zfp281")

ZNF_summary <- ZNF_summary[variable %in% unique(exp_list),]

levels <- ZNF_summary[order(Mean, decreasing=T),]$variable

GEX_ZNF_reshape <- GEX_ZNF_reshape[GEX_ZNF_reshape$variable %in% unique(exp_list),]
GEX_ZNF_reshape$variable <- factor(GEX_ZNF_reshape$variable, levels=levels)

base_plot <- ggplot(GEX_ZNF_reshape[GEX_ZNF_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("ZNF_SCT_GEX-CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="ZNF TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

In the "Ren1_LO" population only.

```R
ZNF_present <- character()
for (gene in ZNF_members) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        ZNF_present <- c(ZNF_present, gene)
    }
}

Ren1_LO_ZNF <- subset(Ren1_LO_GEX, select=c(ZNF_present, "group", "named_cluster"))
Ren1_LO_ZNF_reshape <- reshape2::melt(Ren1_LO_ZNF)

ZNF_summary <- as.data.table(Ren1_LO_ZNF_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))

exp_list <- c(as.character(ZNF_summary[Mean>0 & Median>1,]$variable), "Zfp422", "Zfp467", "Zfp281")

ZNF_summary <- ZNF_summary[variable %in% unique(exp_list),]

levels <- ZNF_summary[order(Mean, decreasing=T),]$variable

Ren1_LO_ZNF_reshape <- Ren1_LO_ZNF_reshape[Ren1_LO_ZNF_reshape$variable %in% unique(exp_list),]
Ren1_LO_ZNF_reshape$variable <- factor(Ren1_LO_ZNF_reshape$variable, levels=levels)

base_plot <- ggplot(Ren1_LO_ZNF_reshape[Ren1_LO_ZNF_reshape$group=="ctrl",], aes(x=variable, y=value)) 
svglite(paste0("ZNF_SCT_CTRL_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot + geom_jitter(width = 0.125, stroke = 0, size=3,
                        aes(col=gene_signature_colors[["Ren1_LO"]],
                        alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="ZNF TF Family") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
```

Plot expression of the subset of overlapping genes

```R
overlap_members <- c("Asap1", "Tbc1d2", "Rsu1", "Tm4sf1", "Btg1",
                     "Ndrg4", "Anxa4", "Il17ra", "Scd1")

overlap_present <- character()
for (gene in c(overlap_members)) {
    if (length(Ren1_LO_GEX[[gene]]) != 0) {
        overlap_present <- c(overlap_present, gene)
    }
}

Ren1_LO_overlap <- subset(Ren1_LO_GEX, select=c(overlap_present, "group", "named_cluster"))
Ren1_LO_overlap_reshape <- reshape2::melt(Ren1_LO_overlap)

overlap_summary <- as.data.table(Ren1_LO_overlap_reshape%>%
group_by(variable)%>% 
summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)))
levels <- overlap_summary[order(Mean, decreasing=T),]$variable
Ren1_LO_overlap_reshape$variable <- factor(Ren1_LO_overlap_reshape$variable, levels=levels)

base_plot1 <- ggplot(Ren1_LO_overlap_reshape[Ren1_LO_overlap_reshape$group == "trt",], aes(x=variable, y=value)) 
base_plot2 <- ggplot(Ren1_LO_overlap_reshape[Ren1_LO_overlap_reshape$group == "ctrl",], aes(x=variable, y=value)) 
svglite(paste0("overlap_trt-ctrl_SCT_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot1 + geom_jitter(width = 0.125, stroke = 0, size=3,
        aes(col=gene_signature_colors[["Ren1_LO"]], alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="H89 Final Gene Set") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
base_plot2 + geom_jitter(width = 0.125, stroke = 0, size=3,
        aes(col=gene_signature_colors[["Ren1_LO"]], alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="DMSO Final Gene Set") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
dev.off()

for (gene in overlap_members) {
    base_plot <- ggplot(Ren1_LO_GEX, aes(x=named_cluster, y=!!rlang::sym(gene)))
    svglite(paste0("final_set_", gene, "_Ren1_LO_SCT_boxplot_",
                   Sys.Date(), ".svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(base_plot + 
        geom_jitter(width = 0.125, stroke = 0, size=3,
                    aes(col=Ren1_LO_GEX$named_cluster,
                        group=Ren1_LO_GEX$named_cluster, alpha=0.1)) +
        geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
        geom_violin(width = 0.50, fill = "transparent",
                    draw_quantiles = 0.5, col='#333333') +
        stat_summary(fun = "mean", geom = "point",
                     shape = 1, size = 2) +
        labs(x="", y="normalized counts", title=gene) +
        scale_color_manual(values=gene_signature_colors) +
        scale_fill_manual(values=gene_signature_colors) +
        custom_theme() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    ylim(0,6)
dev.off()
    dev.off()
} 

for (gene in overlap_members) {
    g1 <- ggstatsplot::ggbetweenstats(
          data  = Ren1_LO_GEX,
          x     = group,
          y     = !!rlang::sym(gene),
          title = gene,
          ggtheme = custom_theme() + theme(panel.grid.major = element_line(
                                            color = "gray",
                                            linewidth = 0.1,
                                            linetype = 2))
        )
    svglite(paste0("final_set_", gene, "_Ren1_LO_SCT_ggbetweenstats_",
            Sys.Date(), ".svg"),
            bg = "transparent", fix_text_size = FALSE, pointsize = 4)
    print(g1)
    dev.off()
}

overlap_present <- character()
for (gene in c(overlap_members)) {
    if (length(GEX_dt[[gene]]) != 0) {
        overlap_present <- c(overlap_present, gene)
    }
}

GEX_overlap <- subset(GEX_dt, select=c(overlap_present, "group", "named_cluster"))
GEX_overlap_reshape <- reshape2::melt(GEX_overlap)

base_plot1 <- ggplot(GEX_overlap_reshape[GEX_overlap_reshape$group == "trt",], aes(x=variable, y=value)) 
base_plot2 <- ggplot(GEX_overlap_reshape[GEX_overlap_reshape$group == "ctrl",], aes(x=variable, y=value)) 
svglite(paste0("overlap-all_trt-ctrl_SCT_boxplot_", Sys.Date(), ".svg"),
               bg = "transparent", fix_text_size = FALSE, pointsize = 4)
base_plot1 + geom_jitter(width = 0.125, stroke = 0, size=3, aes(col=gene_signature_colors[["Ren1_LO"]], alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="H89 Final Gene Set") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    ylim(0,8) +
base_plot2 + geom_jitter(width = 0.125, stroke = 0, size=3, aes(col=gene_signature_colors[["Ren1_LO"]], alpha=0.1)) +
    geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
    geom_violin(width = 0.50, fill = "transparent",
                draw_quantiles = 0.5, col='#333333') +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    labs(x="", y="normalized counts", title="DMSO Final Gene Set") +
    scale_color_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    scale_fill_manual(values=gene_signature_colors[["Ren1_LO"]]) +
    custom_theme() +
    ylim(0,8)
dev.off()
```