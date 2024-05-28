# Plotting qPCR Results 

Set up R environment

```R
library(dplyr)
library(data.table)
library(ggplot2)
library(ggstatsplot)
library(svglite)

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

## A-485

Data from A-485 treatment and washout experiments in As4.1 cells.

```R
A485_qPCR <- fread('As4.1_A485_qPCR_results.csv')

A485_qPCR$condition <- factor(A485_qPCR$condition,
                              levels=c("vehicle",
                                       "0.25uM_trt",
                                       "0.50uM_trt",
                                       "1uM_trt",
                                       "5uM_trt",
                                       "10uM_trt",
                                       "5uM_wash",
                                       "10uM_wash"))

base_plot_A485 <- ggplot(A485_qPCR, aes(x=group, y=Avg_relative_expression)) +
  geom_boxplot() +
  geom_line(aes(col = group, group = group))

svglite("As4.1_A485_qPCR.svg")
base_plot_A485 + geom_jitter(width = 0.125, stroke = 0, size=3, aes(alpha=0.1)) +
  geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
  geom_violin(width = 0.50, fill = "transparent",
              draw_quantiles = 0.5, col='#333333') +
  stat_summary(fun = "mean", geom = "point",
               shape = 1, size = 2) +
  labs(x="", y="Relative Expression to s14", title="As4.1 A485 qPCR") +
  custom_theme() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 0.25,0.50,0.75,1.00,1.25))
dev.off()

A485_statsplot <- ggstatsplot::ggbetweenstats(
  data  = A485_qPCR,
  x     = condition,
  y     = Avg_relative_expression,
  title = "Ren1",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  type = "robust",
  ggtheme = custom_theme() +
    theme(panel.grid.major = element_line(
      color = "gray",
      linewidth = 0.1,
      linetype = 2))
)
svglite("As4.1_A485-all_qPCR.svg")
A485_statsplot
dev.off()

A485_statsplot_5 <- ggstatsplot::ggbetweenstats(
  data  = A485_qPCR[expt == "washout" & group == "5 uM" | group == "vehicle",],
  x     = condition,
  y     = Avg_relative_expression,
  title = "Ren1",
  ggtheme = custom_theme() +
    theme(panel.grid.major = element_line(
      color = "gray",
      linewidth = 0.1,
      linetype = 2))
)
svglite("As4.1_A-485_5uM_qPCR.svg")
A485_statsplot_5 + scale_y_continuous(limits = c(0, 1.25),
                        breaks = c(0, 0.25,0.50,0.75,1.00,1.25))
dev.off()

A485_statsplot_10 <- ggstatsplot::ggbetweenstats(
  data  = A485_qPCR[expt == "washout" & group == "10 uM" | group == "vehicle",],
  x     = condition,
  y     = Avg_relative_expression,
  title = "Ren1",
  ggtheme = custom_theme() +
    theme(panel.grid.major = element_line(
      color = "gray",
      linewidth = 0.1,
      linetype = 2))
)
svglite("As4.1_A-485_10uM_qPCR.svg")
A485_statsplot_10 + scale_y_continuous(limits = c(0, 1.25),
                                    breaks = c(0, 0.25,0.50,0.75,1.00,1.25))
dev.off()
```

## JQ1

Data from JQ1 treatment and washout experiments in As4.1 cells.

```R
JQ1_qPCR <- fread('As4.1_JQ1_qPCR_results.csv')

JQ1_qPCR$condition <- factor(JQ1_qPCR$condition,
                             levels=c("vehicle",
                                      "0.5uM_trt",
                                      "1uM_trt",
                                      "5uM_trt",
                                      "1uM_wash"))
JQ1_statsplot <- ggstatsplot::ggbetweenstats(
  data  = JQ1_qPCR[!grep('^5uM', condition),],
  x     = condition,
  y     = Avg_relative_expression,
  title = "Ren1",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  type = "robust",
  ggtheme = custom_theme() +
    theme(panel.grid.major = element_line(
    color = "gray",
    linewidth = 0.1,
    linetype = 2))
)
svglite("As4.1_JQ1_qPCR.svg")
JQ1_statsplot + scale_y_continuous(limits = c(0, 2),
                        breaks = c(0, 0.25,0.50,0.75,1.00,1.25,1.5,1.75,2.00))
dev.off()
```

## H89

Data from H89 treatment and washout experiments in As4.1 cells.

```R
H89_qPCR <- fread('As4.1_H89_qPCR_results.csv')

H89_qPCR$group <- factor(H89_qPCR$group,
                         levels=c("vehicle_48hr",
                                  "vehicle_96hr",
                                  "vehicle_120hr",
                                  "5uM_48hr_trt",
                                  "5uM_96hr_trt",
                                  "5uM_120hr_trt",
                                  "10uM_48hr_trt",
                                  "10uM_96hr_trt",
                                  "10uM_120hr_trt",
                                  "5uM_96hr_wash",
                                  "5uM_120hr_wash",
                                  "10uM_96hr_wash",
                                  "10uM_120hr_wash"))

base_plot_H89 <- ggplot(H89_qPCR,
                        aes(x=group, y=Avg_relative_expression)) +
  geom_boxplot() +
  geom_line(aes(col = group, group = group))

svglite("As4.1_H89_qPCR.svg")
base_plot_H89 + geom_jitter(width = 0.125, stroke = 0, size=3, aes(alpha=0.1)) +
  geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
  geom_violin(width = 0.50, fill = "transparent",
              draw_quantiles = 0.5, col='#333333') +
  stat_summary(fun = "mean", geom = "point",
               shape = 1, size = 2) +
  labs(x="", y="Relative Expression to s14", title="As4.1 H89 qPCR") +
  custom_theme() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 0.25,0.50,0.75,1.00,1.25))
dev.off()

H89_statsplot <- ggstatsplot::ggbetweenstats(
  data  = H89_qPCR,
  x     = group,
  y     = Avg_relative_expression,
  title = "Ren1",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  type = "robust",
  ggtheme = custom_theme() +
    theme(panel.grid.major = element_line(
      color = "gray",
      linewidth = 0.1,
      linetype = 2))
)
svglite("As4.1_H89-all_qPCR.svg")
H89_statsplot
dev.off()

H89_statsplot_48hr <- ggstatsplot::ggbetweenstats(
  data  = H89_qPCR[grep('48', H89_qPCR$group),],
  x     = group,
  y     = Avg_relative_expression,
  pairwise.comparisons = TRUE,
  pairwise.display = "all",
  type = "robust",
  title = "Ren1",
  ggtheme = custom_theme() +
    theme(panel.grid.major = element_line(
      color = "gray",
      linewidth = 0.1,
      linetype = 2))
)
svglite("As4.1_H89-48hr_qPCR.svg")
H89_statsplot_48hr + scale_y_continuous(limits = c(0, 2.0),
                                        breaks = c(0, 0.25,0.50,0.75,1.00,
                                                   1.25,1.50,1.75,2.00))
dev.off()

H89_96 <- H89_qPCR[grep('96', H89_qPCR$group),]
H89_96$group <- factor(H89_96$group,
                        levels=c("vehicle_96hr",
                                 "5uM_96hr_trt",
                                 "5uM_96hr_wash",
                                 "10uM_96hr_trt",
                                 "10uM_96hr_wash"))  
H89_statsplot_96hr <- ggstatsplot::ggbetweenstats(
  data  = H89_96,
  x     = group,
  y     = Avg_relative_expression,
  pairwise.comparisons = TRUE,
  pairwise.display = "all",
  type = "robust",
  title = "Ren1",
  ggtheme = custom_theme() +
    theme(panel.grid.major = element_line(
      color = "gray",
      linewidth = 0.1,
      linetype = 2))
)
svglite("As4.1_H89-96hr_qPCR.svg")
H89_statsplot_96hr + scale_y_continuous(limits = c(0, 2.0),
                                        breaks = c(0, 0.25,0.50,0.75,1.00,
                                                   1.25,1.50,1.75,2.00))
dev.off()

H89_120 <- H89_qPCR[grep('120', H89_qPCR$group),]
H89_120$group <- factor(H89_120$group,
                        levels=c("vehicle_120hr",
                                 "5uM_120hr_trt",
                                 "5uM_120hr_wash",
                                 "10uM_120hr_trt",
                                 "10uM_120hr_wash"))                    
                    
H89_statsplot_120hr <- ggstatsplot::ggbetweenstats(
  data  = H89_120,
  x     = group,
  y     = Avg_relative_expression,
  pairwise.comparisons = TRUE,
  pairwise.display = "all",
  type = "robust",
  title = "Ren1",
  ggtheme = custom_theme() +
    theme(panel.grid.major = element_line(
      color = "gray",
      linewidth = 0.1,
      linetype = 2))
)
svglite("As4.1_H89-120hr_qPCR.svg")
H89_statsplot_120hr + scale_y_continuous(limits = c(0, 2.0),
                                        breaks = c(0, 0.25,0.50,0.75,1.00,
                                                   1.25,1.50,1.75,2.00))
dev.off()
```

## Combine data to generate a single line graph style plot 

### H89

```R
H89_qPCR <- fread('As4.1_H89_qPCR_results.csv')
H89_qPCR[,inhibitor:="H89"]
H89_qPCR[,time:=c(rep("48", 9), rep("96", 15), rep("120", 15))]
H89_qPCR[,rep:=rep(c("1","2","3"), 13)]
H89_qPCR[,treatment:=c(rep(c("vehicle", "5uM", "10uM"), each=3, times=2),
                       rep(c("5uM", "10uM"), each=3, times=1),
                       rep(c("vehicle", "5uM", "10uM"), each=3, times=1),
                       rep(c("5uM", "10uM"), each=3, times=1))]
H89_10uM <- H89_qPCR[grep('10uM|vehicle', H89_qPCR$group),]
H89_10uM$time <- factor(H89_10uM$time,
                        levels=c("48",
                                 "96",
                                 "120")) 
H89_10uM[,order:=c(rep("1", 3),
                   rep("2", 3),
                   rep("1", 3),
                   rep("5", 3),
                   rep("3", 3),
                   rep("1", 3),
                   rep("6", 3),
                   rep("4", 3))]

H89_line_plot <- ggplot(H89_10uM, aes(x=order, y=Avg_relative_expression)) +
  geom_point(col = 'blue') +
  geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
  geom_violin(width = 0.50, fill = "transparent",
              draw_quantiles = 0.5, col='#333333') +
  stat_summary(fun = "mean", geom = "point",
               shape = 1, size = 2) +
  geom_smooth(aes(x=as.numeric(H89_10uM$order), 
                  y=H89_10uM$Avg_relative_expression), col='blue',
              method=lm, formula = y ~ splines::bs(x, 4), se=FALSE, fullrange=TRUE) +
  labs(x="", y="Relative Expression to s14", title="As4.1 H89 qPCR") +
  scale_y_continuous(limits = c(0, 1.5), breaks = c(0, 0.25,0.50,0.75,1.00,1.25, 1.50)) +
  scale_x_discrete(labels=c("DMSO",
                            "H89 (10uM) - 48hrs",
                            "H89 (10uM) - 96hrs",
                            "H89 (10uM) - 120hrs",
                            "Washout - +96hrs",
                            "Washout - +120hrs")) +
  theme_classic(base_size = 6) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.position = "top",
    legend.justification="center",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )
```

### A485

```R
A485_qPCR <- fread('As4.1_A485_qPCR_results.csv')
A485_10uM <- A485_qPCR[grep('10uM|vehicle', A485_qPCR$treatment),]

A485_line_plot <- ggplot(A485_10uM, aes(x=as.character(order), y=Avg_relative_expression)) +
  geom_point(aes(col = treatment, group = treatment)) +
  geom_boxplot(fill = "transparent", width = 0.5, col='#333333') +
  geom_violin(width = 0.50, fill = "transparent",
              draw_quantiles = 0.5, col='#333333') +
  stat_summary(fun = "mean", geom = "point",
               shape = 1, size = 2) +
  geom_smooth(aes(x=as.numeric(A485_10uM$order),
                  y=A485_10uM$Avg_relative_expression),
              method=lm, formula = y ~ splines::bs(x, 4), se=FALSE, fullrange=TRUE) +
  labs(x="", y="Relative Expression to s14", title="As4.1 A485 qPCR") +
  scale_y_continuous(limits = c(0, 1.5), breaks = c(0, 0.25,0.50,0.75,1.00,1.25, 1.50)) +
  scale_x_discrete(breaks=c(1,2,3,4,5,6),
                   labels=c("DMSO",
                            "A485 (10uM) - 48hrs",
                            "A485 (10uM) - 96hrs",
                            "A485 (10uM) - 120hrs",
                            "Washout - +96hrs",
                            "Washout - +120hrs")) +
  theme_classic(base_size = 6) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.position = "top",
    legend.justification="center",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )
```

### JQ1

```R
JQ1_qPCR <- fread('As4.1_JQ1_qPCR_results.csv')
JQ1_1uM <- JQ1_qPCR[grep('1uM|vehicle', JQ1_qPCR$treatment),]
```

### Plot all treatments together

```R
svglite("As4.1_qPCR_line_plot.svg",
        bg = "transparent", fix_text_size = FALSE, pointsize = 4)
H89_line_plot + geom_point(data=A485_10uM, aes(x=as.character(order), y=Avg_relative_expression,
                                               col = 'red', group = inhibitor)) +
  geom_boxplot(data=A485_10uM, aes(x=as.character(order), y=Avg_relative_expression),
               fill = "transparent", width = 0.5, col='#333333') +
  geom_violin(data=A485_10uM, aes(x=as.character(order), y=Avg_relative_expression),
              width = 0.50, fill = "transparent",
              draw_quantiles = 0.5, col='#333333') +
  stat_summary(fun = "mean", geom = "point",
               shape = 1, size = 2) +
  geom_smooth(aes(x=as.numeric(A485_10uM$order),
                  y=A485_10uM$Avg_relative_expression), col='red',
              method=lm, formula = y ~ splines::bs(x, 4), se=FALSE, fullrange=TRUE) +
  labs(x="", y="Relative Expression to s14", title="As4.1 A485 qPCR") +
  scale_y_continuous(limits = c(0, 1.5), breaks = c(0, 0.25,0.50,0.75,1.00,1.25, 1.50)) +
  geom_point(data=JQ1_1uM, aes(x=as.character(order), y=Avg_relative_expression,
                                   col = 'green', group = inhibitor)) +
  geom_boxplot(data=JQ1_1uM, aes(x=as.character(order), y=Avg_relative_expression),
               fill = "transparent", width = 0.5, col='#333333') +
  geom_violin(data=JQ1_1uM, aes(x=as.character(order), y=Avg_relative_expression),
              width = 0.50, fill = "transparent",
              draw_quantiles = 0.5, col='#333333') +
  stat_summary(fun = "mean", geom = "point",
               shape = 1, size = 2) +
  geom_smooth(aes(x=as.numeric(JQ1_1uM$order),
                  y=JQ1_1uM$Avg_relative_expression), col='green',
              method=lm, formula = y ~ splines::bs(x, 4), se=FALSE, fullrange=TRUE) +
  theme_classic(base_size = 6) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.position = "top",
    legend.justification="center",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )
dev.off()
```

