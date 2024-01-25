### 6. Plotting motif density around ATAC peak summits

```console
cd $PROCESSED/H89/bam_files/FIMO/tomtom/fimo_composites/

for bed in *2M.bed
do
name=$(echo $bed | awk -F"/" '{print $NF}' | awk -F"_2M.bed" '{print $1}')
echo $name
#summing scores of motifs w/in peak
cat $bed | mergeBed -i stdin -c 4 -o sum > ${name}_merged_2M.bed
bedGraphToBigWig ${name}_merged_2M.bed mm10.chrom.sizes ${name}_mm10_instances.bigWig
done
```

```console
#rename 2M files for generating final ATAC dataframe
cp PSWM_family_1_merged_mm10_instances.bigWig main_figure_beds/Fra2_mm10_instances.bigWig
cp PSWM_family_2a_merged_mm10_instances.bigWig main_figure_beds/Klf1_mm10_instances.bigWig
cp PSWM_family_2b_merged_mm10_instances.bigWig main_figure_beds/Sp9_mm10_instances.bigWig
cp PSWM_family_3_merged_mm10_instances.bigWig main_figure_beds/Rorb_mm10_instances.bigWig
cp PSWM_family_4_merged_mm10_instances.bigWig main_figure_beds/Esrrb_mm10_instances.bigWig
cp PSWM_family_5_merged_mm10_instances.bigWig main_figure_beds/ETS-RUNX_mm10_instances.bigWig
cp PSWM_family_6_merged_mm10_instances.bigWig main_figure_beds/Nr5a2_mm10_instances.bigWig
cp PSWM_family_7_merged_mm10_instances.bigWig main_figure_beds/PU.1_mm10_instances.bigWig
cp PSWM_family_8_merged_mm10_instances.bigWig main_figure_beds/RXR_mm10_instances.bigWig
cp PSWM_family_9_merged_mm10_instances.bigWig main_figure_beds/Zfp281_mm10_instances.bigWig
cp PSWM_family_10_merged_mm10_instances.bigWig main_figure_beds/ZNF467_mm10_instances.bigWig
cp PSWM_family_11_merged_mm10_instances.bigWig main_figure_beds/Patz1_mm10_instances.bigWig
cp PSWM_family_12_merged_mm10_instances.bigWig main_figure_beds/Zfp422_mm10_instances.bigWig
cp PSWM_family_13_merged_mm10_instances.bigWig main_figure_beds/Zfp467_mm10_instances.bigWig

cp PSWM_family_ETS_mm10_instances.bigWig main_figure_beds/ETS_mm10_instances.bigWig
cp PSWM_family_ZFPs_mm10_instances.bigWig main_figure_beds/ZFPs_mm10_instances.bigWig
```

```R
library(lattice)
library(bigWig)

dir = '$PROCESSED/H89/bam_files/FIMO/tomtom/fimo_composites/'
setwd(dir)

bedWindow <- function(bed, half_window) {
    bed[,2] = (bed[,2] + bed[,3])/2 - half_window
    bed[,3] = bed[,2] + 2 * half_window
    bed = subset(bed, bed[,2] > 0)
    return(bed)
}

plotFimoLattice <- function(dat, fact = 'Motif', plot_name = 'motif_enrichment_around_summits.pdf',
                            summit = 'Hypersensitivity Summit', class= '',
                            num.m = -200, num.p =90, y.low =0, y.high = 0.2,
                            col_lines = c(rgb(0,0,1,1/2),
                                          rgb(1,0,0,1/2),
                                          rgb(0.1,0.5,0.05,1/2),
                                          rgb(0,0,0,1/2),
                                          rgb(1/2,0,1/2,1/2),
                                          rgb(0,1/2,1/2,1/2),
                                          rgb(1/2,1/2,0,1/2)),
                            fill_poly = c(rgb(0,0,1,1/4),
                                          rgb(1,0,0,1/4),
                                          rgb(0.1,0.5,0.05,1/4),
                                          rgb(0,0,0,1/4),
                                          rgb(1/2,0,1/2,1/4))) {
    
    pdf(plot_name)#, width=6.83, height=3.5)
    print(xyplot(density ~ range|tf, groups = category, data = dat, type = 'l',
                 scales=list(x=list(cex=0.8,relation = "free"),
                             y =list(cex=0.8,axs = 'i',relation = "free")),
                 xlim=c(num.m, num.p),                                       
                 col = col_lines,
                 auto.key = list(points=F, lines=T, cex=0.8, columns = 2),
                 par.settings = list(superpose.symbol = list(pch = c(16),
                                                             col=col_lines,
                                                             cex =0.7),
                                     superpose.line = list(col = col_lines,
                                                           lwd=c(2,2,2,2,2,2),
                                                           lty = c(1,1,1,1,1,1,1,1,1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=0.9, font=1, col='black',font=2),
                 aspect=1.0,
                 between=list(y=0.5, x=0.5),
                 index.cond = list(c(4:6,1:3)),
                 lwd=2,
                 ylab = list(label = "Weighted Motif Density", cex =1,font=2),
                 xlab = list(label = 'Distance from ATAC-seq Peak Summit', cex =1,font=2),
                 strip = function(..., which.panel, bg) {
                     bg.col = 'grey'#c("blue","grey65","red")
                     strip.default(..., which.panel = which.panel,
                                   bg = rep(bg.col, length = which.panel)[which.panel])
                 }
                 ))
    dev.off()
}

fimo_scores_atac <- fimo_scores_all_atac

plot_df = data.frame()
for(i in 1:ncol(fimo_scores_atac)) {
    temp = plot_df_atac[plot_df_atac$peak %in% rownames(fimo_scores_atac[!is.na(fimo_scores_atac[,i]),]),]
    temp$family = colnames(fimo_scores_atac)[i]
    plot_df = rbind(plot_df,temp)
}

plot_df$status = 'Activated'
plot_df[plot_df$supercluster == 'down',]$status = 'Repressed'

all_fimo = data.frame(matrix(ncol = 4, nrow = 0))
colnames(all_fimo) = c('density', 'tf', 'category', 'range')

half_win = 600
file_suffix = '_mm10_instances.bigWig'
dir = '$PROCESSED/H89/bam_files/FIMO/tomtom/fimo_composites/main_figure_beds/'

decreased = plot_df[plot_df$status == 'Repressed',5:7]
decreased[,2] = as.numeric(decreased[,2])
decreased[,3] = as.numeric(decreased[,3])
decreased = bedWindow(decreased, half_win)

increased = plot_df[plot_df$status == 'Activated',5:7]
increased[,2] = as.numeric(increased[,2])
increased[,3] = as.numeric(increased[,3])
increased = bedWindow(increased, half_win)

not_different = read.table('$PROCESSED/H89/bam_files/nondynamic_peaks.bed')
not_different = not_different[not_different$V1 != 'chrM',]
not_different = bedWindow(not_different, half_win)

for(i in 1:ncol(fimo_scores_atac)) {
    factor = colnames(fimo_scores_atac)[i]

    mod_bigWig = paste0(dir, factor, file_suffix)
    factor_name = factor
    print(factor_name)

    loaded_bw = load.bigWig(mod_bigWig)
    
    dec_inten = bed.step.probeQuery.bigWig(loaded_bw, decreased,
                                           gap.value = 0, step = 10,
                                           as.matrix = TRUE)
    dec_query_df = data.frame(cbind(colMeans(dec_inten), factor_name,
                                    'Decreased',
                                    seq(-half_win, (half_win-10), 10)),
                                    stringsAsFactors=F)
    colnames(dec_query_df) = c('density', 'tf', 'category', 'range')
    
    inc_inten = bed.step.probeQuery.bigWig(loaded_bw, increased,
                                           gap.value = 0, step = 10,
                                           as.matrix = TRUE)
    inc_query_df = data.frame(cbind(colMeans(inc_inten), factor_name,
                                    'Increased',
                                    seq(-half_win,(half_win-10), 10)),
                                    stringsAsFactors=F)
    colnames(inc_query_df) = c('density', 'tf', 'category', 'range')
    
    ctrl_inten = bed.step.probeQuery.bigWig(loaded_bw, not_different,
                                            gap.value = 0, step = 10,
                                            as.matrix = TRUE)
    ctrl_query_df = data.frame(cbind(colMeans(ctrl_inten), factor_name,
                                     'Nondynamic',
                                     seq(-half_win, (half_win-10), 10)),
                                     stringsAsFactors=F)
    colnames(ctrl_query_df) = c('density', 'tf', 'category', 'range')
    
    tf_all = rbind(dec_query_df, inc_query_df, ctrl_query_df)

    all_fimo = rbind(all_fimo, tf_all)
}

all_fimo[,1] = as.numeric(all_fimo[,1])
all_fimo[,4] = as.numeric(all_fimo[,4])
```

```
plotFimoLattice(all_fimo[all_fimo$tf %in% c("ETS", "Fra2", "Nr5a2", "PU.1"),],
                num.m = -500, num.p = 500,
                plot_name = 'motif_enrichment_around_summits_decreased-factors.pdf',
                col_lines = c("#FFC20A", "#0C7BDC", "#A2A2A0"))

plotFimoLattice(all_fimo[all_fimo$tf %in% c("Klf1", "Patz1", "Rorb", "Sp9", "ZFPs"),],
                num.m = -500, num.p = 500,
                plot_name = 'motif_enrichment_around_summits_increased-factors.pdf',
                col_lines = c("#FFC20A", "#0C7BDC", "#A2A2A0"))

custom_theme <- function(base_family = "sans", ...){
  theme_classic(base_family = base_family, base_size = 14, ...) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    aspect.ratio = 1,
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )
}

# colors: yellow, blue, gray
# #FFC20A", "#0C7BDC", "#A2A2A0"

# Plot the number of increased or decreased peaks by transcription factor family
svglite("atac_peaks_by_factor.svg",
        bg = "transparent", fix_text_size = FALSE, pointsize = 6)
ggplot(plot_df, aes(x=as.factor(family), fill=as.factor(status) )) +
  geom_bar() +
  scale_fill_manual(name = "Status", labels = c("Increased", "Decreased"),
                    values = c("#FFC20A", "#0C7BDC")) +
  labs(x="", y="Dynamic ATAC peaks") + custom_theme()
dev.off()
```