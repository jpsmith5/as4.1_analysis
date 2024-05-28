## 10. Generate Waddington plot of scRNA-seq cell populations based on gene signatures

Set up R environment.

```R
library(data.table)
library(stringr)
library(rtracklayer)
library(svglite)
library(epiDisplay)
library(magrittr)
library(foreign)
library(psych)
library(tidyverse)
library(rstatix)

custom_theme <- function(base_family = "sans", ...){
  theme_classic(base_family = base_family, base_size = 8, ...) +
    theme(
      axis.line = element_blank(),
      #axis.title.x = element_blank(),
      #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
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
```

Calculate cell numbers and fractions

```R
ctrl_n <- 3114
trt_n  <- 3730

trt_pct <- trt_n/(ctrl_n+trt_n)
ctrl_pct <- ctrl_n/(ctrl_n+trt_n)

all_n <- 3268
remodel_n <- 1290
ren1_lo_n <- 414
khdrbs3_acaa1b_n <- 1183
sec24a_arid3b_n <- 371
gbp7_n <- 318

total_n <- sum(all_n, remodel_n, ren1_lo_n, khdrbs3_acaa1b_n, sec24a_arid3b_n, gbp7_n)
all_pct <- all_n/total_n
remodel_pct <- remodel_n/total_n
ren1_lo_pct <- ren1_lo_n/total_n
khdrbs3_acaa1b_pct <- khdrbs3_acaa1b_n/total_n
sec24a_arid3b_pct <- sec24a_arid3b_n/total_n
gbp7_pct <- gbp7_n/total_n
```

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

## Ctrl samples

```R
all_ctrl_n <- 1511
remodel_ctrl_n <- 561
ren1_lo_ctrl_n <- 90
khdrbs3_acaa1b_ctrl_n <- 588
sec24a_arid3b_ctrl_n <- 192
gbp7_ctrl_n <- 168

total_ctrl_n <- sum(all_ctrl_n, remodel_ctrl_n, ren1_lo_ctrl_n,
  khdrbs3_acaa1b_ctrl_n, sec24a_arid3b_ctrl_n, gbp7_ctrl_n)
all_ctrl_pct <- all_ctrl_n/total_ctrl_n
remodel_ctrl_pct <- remodel_ctrl_n/total_ctrl_n
ren1_lo_ctrl_pct <- ren1_lo_ctrl_n/total_ctrl_n
khdrbs3_acaa1b_ctrl_pct <- khdrbs3_acaa1b_ctrl_n/total_ctrl_n
sec24a_arid3b_ctrl_pct <- sec24a_arid3b_ctrl_n/total_ctrl_n
gbp7_ctrl_pct <- gbp7_ctrl_n/total_ctrl_n
```

## Trt samples

```R
all_trt_n <- 1755
remodel_trt_n <- 728
ren1_lo_trt_n <- 323
khdrbs3_acaa1b_trt_n <- 595
sec24a_arid3b_trt_n <- 177
gbp7_trt_n <- 148

total_trt_n <- sum(all_trt_n, remodel_trt_n, ren1_lo_trt_n,
  khdrbs3_acaa1b_trt_n, sec24a_arid3b_trt_n, gbp7_trt_n)
all_trt_pct <- all_trt_n/total_trt_n
remodel_trt_pct <- remodel_trt_n/total_trt_n
ren1_lo_trt_pct <- ren1_lo_trt_n/total_trt_n
khdrbs3_acaa1b_trt_pct <- khdrbs3_acaa1b_trt_n/total_trt_n
sec24a_arid3b_trt_pct <- sec24a_arid3b_trt_n/total_trt_n
gbp7_trt_pct <- gbp7_trt_n/total_trt_n
```

## Combine

```R
scRNAseq_N <- data.table(signature=c(rep("all+", 1511),
                                     rep("remodel", 561),
                                     rep("Ren1_LO", 90),
                                     rep("Khdrbs3|Acaa1b+", 588),
                                     rep("Sec24a|Arid3b+", 192),
                                     rep("Gbp7+", 168),
                                     rep("all+", 1755),
                                     rep("remodel", 728),
                                     rep("Ren1_LO", 323),
                                     rep("Khdrbs3|Acaa1b+", 595),
                                     rep("Sec24a|Arid3b+", 177),
                                     rep("Gbp7+", 148)),
                         treatment=c(rep("ctrl", total_ctrl_n),
                                     rep("trt", total_trt_n)))

scRNAseq_test <- fisher.test(table(scRNAseq_N$signature, scRNAseq_N$treatment),
    simulate.p.value=TRUE)


scRNA_xtab <- table(scRNAseq_N$signature, scRNAseq_N$treatment)
dimnames(scRNA_xtab) <- list(
  signatures = c("all+", "Gbp7+", "Khdrbs3|Acaa1b+", "remodel", "Ren1_LO", "Sec24a|Arid3b+"),
  treatment = c("ctrl", "trt")
)

fisher_results <- row_wise_fisher_test(scRNA_xtab, p.adjust.method = "fdr", detailed=TRUE)
fwrite(fisher_results, "As4.1_scRNAseq_GeneSignatures_FisherTest_results.csv", sep=",")

signature_palette <- c("#95CF70", "#C90076", "#2ED8C4", 
                       "#873C3E", "#272146", "#2433E1")
names(signature_palette) <- c("all+", "remodel", "Khdrbs3|Acaa1b+", "Ren1_LO", "Sec24a|Arid3b+", "Gbp7+")

# combine plot and statistical test with ggbarstats
library(ggstatsplot)
svglite(paste0("As4.1_scRNAseq_GeneSignature_Fisher_plot.svg"),
        bg = "transparent", fix_text_size = FALSE, pointsize = 4)
ggbarstats(
  scRNAseq_N, signature, treatment, 
  results.subtitle = TRUE,
  type="robust",
  paired=FALSE,
  ggtheme=custom_theme()
) + scale_fill_manual(values=signature_palette)
dev.off()

ren1_vs_others <- data.table(signature=c(rep("other", 1511),
                                         rep("other", 561),
                                         rep("Ren1_LO", 90),
                                         rep("other", 588),
                                         rep("other", 192),
                                         rep("other", 168),
                                         rep("other", 1755),
                                         rep("other", 728),
                                         rep("Ren1_LO", 323),
                                         rep("other", 595),
                                         rep("other", 177),
                                         rep("other", 148)),
                             treatment=c(rep("ctrl", total_ctrl_n),
                                         rep("trt", total_trt_n)))
ren1lo_test <- fisher.test(table(ren1_vs_others$signature, ren1_vs_others$treatment))
svglite(paste0("As4.1_scRNAseq_GeneSignature_Fisher-Ren1LO_plot.svg"),
        bg = "transparent", fix_text_size = FALSE, pointsize = 4)
ggbarstats(
  ren1_vs_others, treatment, signature, 
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(ren1lo_test$p.value < 0.001, "< 0.001", round(ren1lo_test$p.value, 3))
  )
)
dev.off()

svglite(paste0("As4.1_scRNAseq_waddington_plot.svg"),
        bg = "transparent", fix_text_size = FALSE, pointsize = 4)
waddingtonplot::waddingtonPlot(branches=c(1,2,12),
                               horizontal_distribution=list(c(1),
                                 c(ctrl_pct,trt_pct),
                                 c(all_ctrl_pct, remodel_ctrl_pct,
                                   ren1_lo_ctrl_pct, khdrbs3_acaa1b_ctrl_pct,
                                   sec24a_arid3b_ctrl_pct, gbp7_ctrl_pct,
                                   all_trt_pct, remodel_trt_pct,
                                   ren1_lo_trt_pct, khdrbs3_acaa1b_trt_pct,
                                   sec24a_arid3b_trt_pct, gbp7_trt_pct)),
                               ridge.count = 75,
                               line.type = 3,
                               do.return = T) + custom_theme()
dev.off()
```
