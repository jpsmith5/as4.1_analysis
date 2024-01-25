### 6. Create heatmap for given feature/peak ranges

#### 6a. H3K27Ac

```R
project_name <- "H3K27Ac"

H3K27Ac_path <- "$PROCESSED/nf-core/chipseq/H89/bwa/mergedLibrary/macs2/broadPeak/"
H3K27Ac_files <- dir(H3K27Ac_path, "broadPeak")
H3K27Ac_data  <- sapply(file.path(H3K27Ac_path, H3K27Ac_files),
                        toGRanges, format="broadPeak")
names(H3K27Ac_data) <- gsub(".broadPeak", "", H3K27Ac_files)

ol       <- ChIPpeakAnno::findOverlapsOfPeaks(H3K27Ac_data)
features <- ol$peaklist[[length(ol$peaklist)]]
wid      <- width(features)
feature.recentered        <- feature.center <- features
start(feature.center)     <- start(features) + floor(wid/2)
width(feature.center)     <- 1
start(feature.recentered) <- start(feature.center) - 2000
end(feature.recentered)   <- end(feature.center) + 2000

bigWig_path <- "$PROCESSED/nf-core/chipseq/H89/bwa/mergedLibrary/bigwig/"
H3K27Ac_bw_files <- dir(bigWig_path, "bigWig")
H3K27Ac_bw_files <- H3K27Ac_bw_files[-5] # drop pooled input

cvglists <- sapply(file.path(bigWig_path, H3K27Ac_bw_files), import,
                   format="BigWig",
                   which=feature.recentered,
                   as="RleList")
names(cvglists) <- gsub(".bigWig", "", H3K27Ac_bw_files)

sig <- featureAlignedSignal(cvglists, feature.center,
                            upstream=2000, downstream=2000)
svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_signal-heatmap.svg"))
heatmap <- featureAlignedHeatmap(sig, feature.center,
                                 upstream=2000, downstream=2000,
                                 upper.extreme=c(3,.5,4))
dev.off()
svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_aligned-signal.svg"))
featureAlignedDistribution(sig, feature.center,
                           upstream=2000, downstream=2000,
                           type="l")
dev.off()
```

#### 6b. P300
```R
project_name <- "P300"

P300_path <- "$PROCESSED/nf-core/chipseq/P300/bwa/mergedLibrary/macs2/broadPeak/"
P300_files <- dir(P300_path, "broadPeak")
P300_data  <- sapply(file.path(P300_path, P300_files),
                        toGRanges, format="broadPeak")
names(P300_data) <- gsub(".broadPeak", "", P300_files)

ol       <- ChIPpeakAnno::findOverlapsOfPeaks(P300_data)
features <- ol$peaklist[[length(ol$peaklist)]]
wid      <- width(features)
feature.recentered        <- feature.center <- features
start(feature.center)     <- start(features) + floor(wid/2)
width(feature.center)     <- 1
start(feature.recentered) <- start(feature.center) - 2000
end(feature.recentered)   <- end(feature.center) + 2000

bigWig_path <- "$PROCESSED/nf-core/chipseq/P300/bwa/mergedLibrary/bigwig/"
P300_bw_files <- dir(bigWig_path, "bigWig")
P300_bw_files <- P300_bw_files[-1] # drop pooled input

cvglists <- sapply(file.path(bigWig_path, P300_bw_files), import,
                   format="BigWig",
                   which=feature.recentered,
                   as="RleList")
names(cvglists) <- gsub(".bigWig", "", P300_bw_files)

sig <- featureAlignedSignal(cvglists, feature.center,
                            upstream=2000, downstream=2000)
svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_signal-heatmap.svg"))
heatmap <- featureAlignedHeatmap(sig, feature.center,
                                 upstream=2000, downstream=2000,
                                 upper.extreme=c(3,.5,4))
dev.off()
svglite(paste0(project_name, "_DB/", project_name,
               "_H89-v-DMSO_aligned-signal.svg"))
featureAlignedDistribution(sig, feature.center,
                           upstream=2000, downstream=2000,
                           type="l")
dev.off()
```