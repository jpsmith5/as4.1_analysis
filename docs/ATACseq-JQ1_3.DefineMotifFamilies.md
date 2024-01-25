## 3. Motif analysis

### 3a. Find motifs in the dynamic peak clusters.

 - cluster_bed_cluster1.bed
 - cluster_bed_cluster2.bed
 
```console
intersectBed -v -a all_peaks.bed -b cluster_bed_cluster1.bed > nondynamic_peaks_cluster1.bed
intersectBed -v -a all_peaks.bed -b cluster_bed_cluster2.bed > nondynamic_peaks_cluster2.bed

mkdir FIMO/
cd FIMO/

## cluster1
slopBed -i ../cluster_bed_cluster1.bed \
      -g mm10.chrom.sizes -b -50 | \
        fastaFromBed -fi mm10.fa -bed stdin \
              -fo cluster1.fasta
cd ../
meme -p 16 -oc FIMO/cluster1_meme -nmotifs 15 -objfun classic -evt 0.01 -searchsize 0 -minw 6 \
         -maxw 18 -revcomp -dna -markov_order 3 -maxsize 100000000 \
         FIMO/cluster1.fasta

## cluster2
cd FIMO/
slopBed -i ../cluster_bed_cluster2.bed \
      -g mm10.chrom.sizes -b -50 | \
        fastaFromBed -fi mm10.fa -bed stdin \
              -fo cluster2.fasta
cd ../
meme -p 1 -oc FIMO/cluster2_meme -nmotifs 15 -objfun classic -evt 0.01 -searchsize 0 -minw 6 \
         -maxw 18 -revcomp -dna -markov_order 3 -maxsize 100000000 \
         FIMO/cluster2.fasta
```

### 3b. Generate combined motif database

Use the one already generated for the H89 analysis.

```console
cd FIMO/
mkdir databases/
cd databases/

cp $PROCESSED/H89/bam_files/FIMO/databases/motif_database_combined.txt ./

cd ../
```

### 3c. Run TomTom on Meme output

```console
mkdir tomtom/
cd tomtom/

echo 'Running TOMTOM'

tomtom -no-ssc -oc cluster1.tomtom_output -verbosity 1 -min-overlap 5 -dist ed -evalue -thresh 0.05 ../cluster1_meme/meme.txt ../databases/motif_database_combined.txt

tomtom -no-ssc -oc cluster2.tomtom_output -verbosity 1 -min-overlap 5 -dist ed -evalue -thresh 0.05 ../cluster2_meme/meme.txt ../databases/motif_database_combined.txt

tomtom -no-ssc -oc cluster1.tomtom_output -verbosity 1 -min-overlap 5 -dist ed -text -evalue -thresh 0.05 ../cluster1_meme/meme.txt ../databases/motif_database_combined.txt > cluster1.tomtom_output/tomtom.txt

tomtom -no-ssc -oc cluster2.tomtom_output -verbosity 1 -min-overlap 5 -dist ed -text -evalue -thresh 0.05 ../cluster2_meme/meme.txt ../databases/motif_database_combined.txt > cluster2.tomtom_output/tomtom.txt
```

 - extract motif names from tomtom output
 - extract individual meme files from combined database

```console
mkdir individual_memes
cd individual_memes/

cp $PROCESSED/H89/bam_files/FIMO/tomtom/individual_memes/MEME_individual_from_db.py .
python2.7 MEME_individual_from_db.py -i ../../databases/motif_database_combined.txt

for file in *meme.txt 
do
    name=$(echo $file | awk -F"motif_database_combined.txt_" '{print $1}')
    mv $file ${name}meme.txt
    
done

#tomtom all query factors against all others
cd ../
mkdir tomtom_all_query_factors
cd tomtom_all_query_factors/
```

`motif_names.txt` is a single column list of the motifs discovered from the above `tomtom` command (lines 633-639). Make one for cluster1 and cluster2.

Copy the ID column from tomtom.tsv: `cluster1_motif_names.txt`, `cluster2_motif_names.txt`

| cluster1                                 | cluster2                                 |
|------------------------------------------|------------------------------------------|
| Klf4_cisbp                               | AP-1_homer                               |
| Klf1_cisbp                               | (Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
| KLF5_homer                               | BATF_homer                               |
| (Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp  | Fos_homer                                |
| KLF1_homer                               | Atf3_homer                               |
| (Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp   | Fra2_homer                               |
| KLF3_homer                               | Fra1_homer                               |
| (Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp   | Jun-AP1_homer                            |
| EKLF_homer                               | Bach1_homer                              |
| KLF14_homer                              | (Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp  |
| KLF6_homer                               | (Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp |
| Klf7_cisbp                               | JunB_homer                               |
| (Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp   | NF-E2_homer                              |
| Klf12_cisbp                              | Bach2_homer                              |
| Klf7_primary_uniprobe                    | Nrf2_homer                               |
| (Klf14)_(Homo_sapiens)_(DBD_0.97)_cisbp  | (Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp   |
| (Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp    | NFE2L2_homer                             |
| (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp    | Fosl2_homer                              |
| Klf4_homer                               | Nfe2l2_cisbp                             |
| Klf8_cisbp                               | MafK_homer                               |
| Sp2_homer                                | (Maff)_(Homo_sapiens)_(DBD_1.00)_cisbp   |
| (Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp  | MafB_homer                               |
| (Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp  | Jund_cisbp                               |
| Klf9_homer                               | ERE_homer                                |
| (Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp    | Nr2c1_cisbp                              |
| (Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp   | Rorb_cisbp                               |
| (Klf15)_(Homo_sapiens)_(DBD_1.00)_cisbp  | Rxrb_cisbp                               |
| KLF10_homer                              | Nr2f2_cisbp                              |
| TEAD2_homer                              | Rarb_cisbp                               |
| (Tead1)_(Homo_sapiens)_(DBD_1.00)_cisbp  | Rara_cisbp                               |
| TEAD1_homer                              | Nr2f6_cisbp                              |
| (Tead4)_(Homo_sapiens)_(DBD_0.94)_cisbp  | Rxrg_cisbp                               |
| TEAD4_homer                              | Esr1_cisbp                               |
| TEAD_homer                               | RUNX1_homer                              |
| (Tead2)_(Homo_sapiens)_(DBD_0.97)_cisbp  | RUNX_homer                               |
| TEAD3_homer                              | RUNX-AML_homer                           |
| BATF_homer                               | RUNX2_homer                              |
| Atf3_homer                               | Runx1_cisbp                              |
| Jun-AP1_homer                            |                                          |
| AP-1_homer                               |                                          |
| Fra2_homer                               |                                          |
| Fos_homer                                |                                          |
| Fra1_homer                               |                                          |
| (Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp    |                                          |
| Bach1_homer                              |                                          |
| JunB_homer                               |                                          |
| Bach2_homer                              |                                          |
| (Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp |                                          |
| Nrf2_homer                               |                                          |
| (Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp  |                                          |
| NF-E2_homer                              |                                          |
| Fosl2_homer                              |                                          |
| (Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp   |                                          |
| MafK_homer                               |                                          |
| NFE2L2_homer                             |                                          |
| MafB_homer                               |                                          |
| Nfe2l2_cisbp                             |                                          |
| (Batf)_(Homo_sapiens)_(DBD_0.98)_cisbp   |                                          |

```console
dos2unix ../cluster1.tomtom_output/cluster1_motif_names.txt
dos2unix ../cluster2.tomtom_output/cluster2_motif_names.txt

cat ../cluster1.tomtom_output/cluster1_motif_names.txt | while read factor
do
    echo $factor
    cp ../individual_memes/${factor}_meme.txt $PWD
done

cat *meme.txt > ../cluster1_query_factors_meme.txt
rm *meme.txt

cat ../cluster2.tomtom_output/cluster2_motif_names.txt | while read factor
do
    echo $factor
    cp ../individual_memes/${factor}_meme.txt $PWD
done

cat *meme.txt > ../cluster2_query_factors_meme.txt

cd ../
cat cluster1_query_factors_meme.txt cluster2_query_factors_meme.txt > all_query_factors_meme.txt
cd tomtom_all_query_factors
```

### 3d. Define motif families

On any of the motifs with a parenthesis, may need to run `tomtom` manually. For example:

| cluster1                                 | cluster2                                 |
|------------------------------------------|------------------------------------------|
| (Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp  | (Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
| (Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp   | (Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp  |
| (Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp   | (Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp |
| (Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp   | (Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp   |
| (Klf14)_(Homo_sapiens)_(DBD_0.97)_cisbp  | (Maff)_(Homo_sapiens)_(DBD_1.00)_cisbp   |
| (Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp    |                                          |
| (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp    |                                          |
| (Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp  |                                          |
| (Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp  |                                          |
| (Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp    |                                          |
| (Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp   |                                          |
| (Klf15)_(Homo_sapiens)_(DBD_1.00)_cisbp  |                                          |
| (Tead1)_(Homo_sapiens)_(DBD_1.00)_cisbp  |                                          |
| (Tead4)_(Homo_sapiens)_(DBD_0.94)_cisbp  |                                          |
| (Tead2)_(Homo_sapiens)_(DBD_0.97)_cisbp  |                                          |
| (Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp    |                                          |
| (Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp |                                          |
| (Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp  |                                          |
| (Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp   |                                          |
| (Batf)_(Homo_sapiens)_(DBD_0.98)_cisbp   |                                          |

```console
mkdir tmp/
mv '(Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp_meme.txt' tmp/
mv '(Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp_meme.txt' tmp/
mv '(Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp_meme.txt' tmp/
mv '(Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Klf14)_(Homo_sapiens)_(DBD_0.97)_cisbp_meme.txt' tmp/
mv '(Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp_meme.txt' tmp/
mv '(Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp_meme.txt' tmp/
mv '(Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Klf15)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Tead1)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Tead4)_(Homo_sapiens)_(DBD_0.94)_cisbp_meme.txt' tmp/
mv '(Tead2)_(Homo_sapiens)_(DBD_0.97)_cisbp_meme.txt' tmp/
mv '(Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp_meme.txt' tmp/
mv '(Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp_meme.txt' tmp/
mv '(Batf)_(Homo_sapiens)_(DBD_0.98)_cisbp_meme.txt' tmp/

mv '(Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp_meme.txt' tmp/
mv '(Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/
mv '(Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp_meme.txt' tmp/
mv '(Maff)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' tmp/

for meme in *meme.txt
do
    name=$(echo $meme | awk -F".txt_" '{print $NF}' | awk -F"_meme.txt" '{print $1}')
    echo $name
    tomtom -no-ssc -oc $name.tomtom_output -verbosity 1  -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 ${meme} ../all_query_factors_meme.txt
    cd $name.tomtom_output
    cut -f1,2,5 tomtom.tsv | tail -n +2 | sed '$d' | sed '$d' | sed '$d' | sed '$d' >> ../3_col_combined_motif_db_pre.txt
    cd ..
done

# Manually run tomtom on odd name files. For example:
tomtom -no-ssc -oc '(Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp.tomtom_output' -verbosity 1  -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 'tmp/(Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp_meme.txt' ../all_query_factors_meme.txt
cd '(Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp.tomtom_output' && cut -f1,2,5 tomtom.tsv | tail -n +2 | sed '$d' | sed '$d' | sed '$d' | sed '$d' >> ../3_col_combined_motif_db_pre.txt && cd ../

grep -v '#' 3_col_combined_motif_db_pre.txt > 3_col_combined_motif_db.txt
rm 3_col_combined_motif_db_pre.txt
cp 3_col_combined_motif_db.txt ..
cd ..
```

visualize pswms

```R
library(igraph)
library(dichromat)

set.seed(99)
setwd('$PROCESSED/jq1/bam_files/FIMO/tomtom/')

threecol=read.table("3_col_combined_motif_db.txt",header=F,stringsAsFactors = F,sep='\t')
colnames(threecol)=c('from','to','e_value')
threecol$weight=abs(log(threecol$e_value))

#create the graph variable
g=graph.data.frame(threecol,directed=F)
g=simplify(g)

cluster=clusters(g)

for(i in 1:length(groups(cluster))) {
write.table(groups(cluster)[i],file=paste0('PSWM_family_',i,'.txt'),col.names = F, row.names = F, quote = F, sep = '\t')
}

l=layout.fruchterman.reingold(g)
l=layout.norm(l,-1,1,-1,1)

colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue","purple"))
#pick a distinct color from the palette for each disease and save the list of colors as a vector
mycol = colfunc(length(groups(cluster)))

pdf(file='families.pdf',width=10,height=10)
plot(g,layout=l,rescale=F,vertex.label.cex=.5,xlim=range(l[,1]),  ylim=range(l[,2]),
edge.width=E(g)$weight/20,vertex.size=degree(g,mode='out')/5,
edge.curved=T,vertex.label=NA,vertex.color=mycol[cluster$membership],
margin=0,asp=0)
dev.off()
```

Split the SP/KLF family into individual families.

```console
cd /scratch/jps3dp/processed/pepatac/gomez/robert/H89/bam_files/FIMO/tomtom/
grep -e '[Kk][Ll][Ff]' PSWM_family_2.txt > PSWM_family_2a.txt
grep -e '[Ss][Pp]' PSWM_family_2.txt > PSWM_family_2b.txt
```

| Family1                                  | Family2                                 | Family2a                                | Family2b                              | Family3     | Family4        | Family5                                 |
|------------------------------------------|-----------------------------------------|-----------------------------------------|---------------------------------------|-------------|----------------|-----------------------------------------|
| AP-1_homer                               | EKLF_homer                              | EKLF_homer                              | Sp2_homer                             | ERE_homer   | Runx1_cisbp    | TEAD1_homer                             |
| Atf3_homer                               | KLF10_homer                             | KLF10_homer                             | (Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp | Esr1_cisbp  | RUNX1_homer    | (Tead1)_(Homo_sapiens)_(DBD_1.00)_cisbp |
| Bach1_homer                              | (Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp | (Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp | (Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp | Nr2c1_cisbp | RUNX2_homer    | TEAD2_homer                             |
| Bach2_homer                              | Klf12_cisbp                             | Klf12_cisbp                             | (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp | Nr2f2_cisbp | RUNX-AML_homer | (Tead2)_(Homo_sapiens)_(DBD_0.97)_cisbp |
| (Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp  | KLF14_homer                             | KLF14_homer                             |                                       | Nr2f6_cisbp | RUNX_homer     | TEAD3_homer                             |
| BATF_homer                               | (Klf14)_(Homo_sapiens)_(DBD_0.97)_cisbp | (Klf14)_(Homo_sapiens)_(DBD_0.97)_cisbp |                                       | Rara_cisbp  |                | TEAD4_homer                             |
| (Batf)_(Homo_sapiens)_(DBD_0.98)_cisbp   | (Klf15)_(Homo_sapiens)_(DBD_1.00)_cisbp | (Klf15)_(Homo_sapiens)_(DBD_1.00)_cisbp |                                       | Rarb_cisbp  |                | (Tead4)_(Homo_sapiens)_(DBD_0.94)_cisbp |
| Fos_homer                                | (Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp | (Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp |                                       | Rorb_cisbp  |                | TEAD_homer                              |
| (Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp    | Klf1_cisbp                              | Klf1_cisbp                              |                                       | Rxrb_cisbp  |                |                                         |
| Fosl2_homer                              | KLF1_homer                              | KLF1_homer                              |                                       | Rxrg_cisbp  |                |                                         |
| Fra1_homer                               | (Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp  | (Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp  |                                       |             |                |                                         |
| Fra2_homer                               | KLF3_homer                              | KLF3_homer                              |                                       |             |                |                                         |
| Jun-AP1_homer                            | (Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp  | (Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp  |                                       |             |                |                                         |
| JunB_homer                               | Klf4_cisbp                              | Klf4_cisbp                              |                                       |             |                |                                         |
| Jund_cisbp                               | Klf4_homer                              | Klf4_homer                              |                                       |             |                |                                         |
| MafB_homer                               | KLF5_homer                              | KLF5_homer                              |                                       |             |                |                                         |
| (Maff)_(Homo_sapiens)_(DBD_1.00)_cisbp   | (Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp  | (Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp  |                                       |             |                |                                         |
| MafK_homer                               | KLF6_homer                              | KLF6_homer                              |                                       |             |                |                                         |
| NF-E2_homer                              | (Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp  | (Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp  |                                       |             |                |                                         |
| (Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp   | Klf7_cisbp                              | Klf7_cisbp                              |                                       |             |                |                                         |
| (Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp | Klf7_primary_uniprobe                   | Klf7_primary_uniprobe                   |                                       |             |                |                                         |
| Nfe2l2_cisbp                             | Klf8_cisbp                              | Klf8_cisbp                              |                                       |             |                |                                         |
| NFE2L2_homer                             | Klf9_homer                              | Klf9_homer                              |                                       |             |                |                                         |
| Nrf2_homer                               | (Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp |                                         |                                       |             |                |                                         |
|                                          | Sp2_homer                               |                                         |                                       |             |                |                                         |
|                                          | (Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp   |                                         |                                       |             |                |                                         |
|                                          | (Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp   |                                         |                                       |             |                |                                         |
|                                          | (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp   |                                         |                                       |             |                |                                         |

```console
cp $PROCESSED/H89/bam_files/FIMO/tomtom/tomtom_output_to_composite.py ./
cp $PROCESSED/H89/bam_files/FIMO/tomtom/meme_header.txt ./
cp $PROCESSED/H89/bam_files/FIMO/tomtom/generate_composite_motif.R ./

mkdir composite_motifs
cd composite_motifs/

mkdir PSWM_family_1/
mkdir PSWM_family_2/
mkdir PSWM_family_2a/
mkdir PSWM_family_2b/
mkdir PSWM_family_3/
mkdir PSWM_family_4/
mkdir PSWM_family_5/
```

#### PSWM Family 1

```console
cd PSWM_family_1/
cat ../../PSWM_family_1.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_1_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc AP-1_homer.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 AP-1_homer_meme.txt PSWM_family_1_factors_meme.txt

max_motif=$(wc -l < AP-1_homer.tomtom_output/tomtom.tsv)
final_query="AP-1_homer"

# Check the rest
echo FINAL_QUERY IS $final_query
echo MAX_MOTIF IS $max_motif
# FINAL_QUERY IS Bach2_homer

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_1_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_1_composite_index.txt
```

Now, for each additional family, do the same process.

#### PSWM Family 2

```console
cd ../../PSWM_family_2/
cat ../../PSWM_family_2.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_2_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'EKLF_homer.tomtom_output' -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 'EKLF_homer_meme.txt' PSWM_family_2_factors_meme.txt && max_motif=$(wc -l < 'EKLF_homer.tomtom_output/tomtom.tsv') && final_query="EKLF_homer"

# Check the rest
echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS Klf1_cisbp

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_2_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_2_composite_index.txt
``` 

##### PSWM Family 2a (KLF)

```console
cd ../../PSWM_family_2a/
cat ../../PSWM_family_2a.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_2a_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'EKLF_homer.tomtom_output' -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 'EKLF_homer_meme.txt' PSWM_family_2a_factors_meme.txt && max_motif=$(wc -l < 'EKLF_homer.tomtom_output/tomtom.tsv') && final_query="EKLF_homer"

# Check the rest
echo FINAL_QUERY IS $final_query
echo MAX_MOTIF IS $max_motif
# FINAL_QUERY IS Klf1_cisbp (max_motif: 20)

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_2a_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_2a_composite_index.txt
#cp ../PSWM_family_2a_composite_values.txt ../composite.values.txt
#cp ../PSWM_family_2a_composite_index.txt ../composite.index.txt
``` 

##### PSWM Family 2b (SP)

```console
cd ../../PSWM_family_2b/
cat ../../PSWM_family_2b.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_2b_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'Sp2_homer.tomtom_output' -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 'Sp2_homer_meme.txt' PSWM_family_2b_factors_meme.txt && max_motif=$(wc -l < 'Sp2_homer.tomtom_output/tomtom.tsv') && final_query="Sp2_homer"

# Check the rest
echo FINAL_QUERY IS $final_query
echo MAX_MOTIF IS $max_motif
# FINAL_QUERY IS (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp (max_motif: 8)

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_2b_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_2b_composite_index.txt
#cp ../PSWM_family_2b_composite_values.txt ../composite.values.txt
#cp ../PSWM_family_2b_composite_index.txt ../composite.index.txt
``` 

#### PSWM Family 3

```console
cd ../../PSWM_family_3/
cat ../../PSWM_family_3.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_3_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'ERE_homer.tomtom_output' -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 'ERE_homer_meme.txt' PSWM_family_3_factors_meme.txt && max_motif=$(wc -l < 'ERE_homer.tomtom_output/tomtom.tsv') && final_query="ERE_homer"

# Check the rest
echo FINAL_QUERY IS $final_query
echo MAX_MOTIF IS $max_motif
# FINAL_QUERY IS Rorb_cisbp (max motif: 15)

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_3_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_3_composite_index.txt
``` 

#### PSWM Family 4

```console
cd ../../PSWM_family_4/
cat ../../PSWM_family_4.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_4_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'Runx1_cisbp.tomtom_output' -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 'Runx1_cisbp_meme.txt' PSWM_family_4_factors_meme.txt && max_motif=$(wc -l < 'Runx1_cisbp.tomtom_output/tomtom.tsv') && final_query="Runx1_cisbp"

echo FINAL_QUERY IS $final_query
echo MAX_MOTIF IS $max_motif
# FINAL_QUERY IS RUNX_homer (max motif: 10)

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_4_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_4_composite_index.txt
``` 

#### PSWM Family 5

```console
cd ../../PSWM_family_5/
cat ../../PSWM_family_5.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_5_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'TEAD1_homer.tomtom_output' -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 'TEAD1_homer_meme.txt' PSWM_family_5_factors_meme.txt && max_motif=$(wc -l < 'TEAD1_homer.tomtom_output/tomtom.tsv') && final_query="TEAD1_homer"

echo FINAL_QUERY IS $final_query
echo MAX_MOTIF IS $max_motif
# FINAL_QUERY IS TEAD4_homer (max motif: 13)

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_5_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_5_composite_index.txt
``` 

 - Must update the directory in `generate_composite_motif.R`

```
cd $PROCESSED/jq1/bam_files/FIMO/tomtom/composite_motifs/

cp PSWM_family_1/PSWM_family_1_composite_index.txt PSWM_family_1/composite.index.txt
cp PSWM_family_1/PSWM_family_1_composite_values.txt PSWM_family_1/composite.values.txt        
cp PSWM_family_2a/PSWM_family_2a_composite_index.txt PSWM_family_2a/composite.index.txt
cp PSWM_family_2a/PSWM_family_2a_composite_values.txt PSWM_family_2a/composite.values.txt
cp PSWM_family_2b/PSWM_family_2b_composite_index.txt PSWM_family_2b/composite.index.txt
cp PSWM_family_2b/PSWM_family_2b_composite_values.txt PSWM_family_2b/composite.values.txt
cp PSWM_family_3/PSWM_family_3_composite_values.txt PSWM_family_3/composite.values.txt
cp PSWM_family_3/PSWM_family_3_composite_index.txt PSWM_family_3/composite.index.txt
cp PSWM_family_4/PSWM_family_4_composite_values.txt PSWM_family_4/composite.values.txt
cp PSWM_family_4/PSWM_family_4_composite_index.txt PSWM_family_4/composite.index.txt
cp PSWM_family_5/PSWM_family_5_composite_values.txt PSWM_family_5/composite.values.txt
cp PSWM_family_5/PSWM_family_5_composite_index.txt PSWM_family_5/composite.index.txt

for family in PSWM*
do
    cd $family
    num=$(ls *txt | wc -l)

    if [[ $num -ge 2 ]]
    then
    
    #generate composite PSWM
    Rscript ../../generate_composite_motif.R $family
    cat ../../meme_header.txt ${family}_composite_PSWM.txt > ${family}_meme.txt	
    else
    line=`grep MOTIF *meme.txt`
    cp *meme.txt ${family}_meme.txt
    sed -i "s;${line};MOTIF   Composite;g" ${family}_meme.txt
    fi
    cd ..
done
```

| PSWM_family | composite_motif                          |
|-------------|------------------------------------------|
| 1           | Bach2_homer                              |
| 2a          | Klf1_cisbp                               |
| 2b          | (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
| 3           | Rorb_cisbp                               |
| 4           | RUNX_homer                               |
| 5           | TEAD4_homer                              |

### 3e. Run FIMO against composite PSWM and take top 2 million hits

```console
cd $PROCESSED/jq1/bam_files/FIMO/tomtom
mkdir fimo_composites/
cd fimo_composites/
```

```console
wget  -O mm10.fa http://awspds.refgenie.databio.org/refgenomes.databio.org/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/fasta__default/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1.fa
```

Family 1 (Bach2_homer)

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_1/PSWM_family_1_meme.txt mm10.fa > PSWM_family_1_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_1_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_1_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_1_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_1_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_1_composite_fimo.txt
```

Family 2a (Klf1_cisbp)

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_2a/PSWM_family_2a_meme.txt mm10.fa > PSWM_family_2a_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_2a_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_2a_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_2a_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_2a_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_2a_composite_fimo.txt
```

Family 2b (Sp9)

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_2b/PSWM_family_2b_meme.txt mm10.fa > PSWM_family_2b_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_2b_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_2b_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_2b_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_2b_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_2b_composite_fimo.txt
```

Family 3 (Rorb_cisbp)

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_3/PSWM_family_3_meme.txt mm10.fa > PSWM_family_3_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_3_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_3_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_3_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_3_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_3_composite_fimo.txt
```

Family 4 (RUNX_homer)

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_4/PSWM_family_4_meme.txt mm10.fa > PSWM_family_4_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_4_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_4_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_4_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_4_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_4_composite_fimo.txt
```

Family 5 (TEAD4_homer)

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_5/PSWM_family_5_meme.txt mm10.fa > PSWM_family_5_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_5_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_5_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_5_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_5_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_5_composite_fimo.txt
```

```console
for i in *_2M.txt
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_2M.txt" '{print $1}')
    echo $name
    intersectBed -loj -a ../../../dynamic_peaks.bed -b $i > ${name}_fimo.bed
    intersectBed -loj -a ../../../nondynamic_peaks.bed -b $i > ${name}_fimo_nondyn.bed
    intersectBed -loj -a ../../../all_peaks.bed -b $i > ${name}_fimo_all.bed
    cat $i | cut -f1-3,5 | sort -k1,1 -k2,2n > ${name}_2M.bed 
done

mkdir main_figure_beds/
```

| PSWM_family | composite_motif                          |
|-------------|------------------------------------------|
| 1           | Bach2_homer                              |
| 2a          | Klf1_cisbp                               |
| 2b          | (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
| 3           | Rorb_cisbp                               |
| 4           | RUNX_homer                               |
| 5           | TEAD4_homer                              |

Copy BED files into figure folder

```console
cp PSWM_family_1_fimo.bed main_figure_beds/Bach2_fimo.bed
cp PSWM_family_2a_fimo.bed main_figure_beds/Klf1_fimo.bed
cp PSWM_family_2b_fimo.bed main_figure_beds/Sp9_fimo.bed
cp PSWM_family_3_fimo.bed main_figure_beds/Rorb_fimo.bed
cp PSWM_family_4_fimo.bed main_figure_beds/RUNX_fimo.bed
cp PSWM_family_5_fimo.bed main_figure_beds/TEAD4_fimo.bed

cp PSWM_family_1_2M.bed main_figure_beds/Bach2_2M.bed
cp PSWM_family_2a_2M.bed main_figure_beds/Klf1_2M.bed
cp PSWM_family_2b_2M.bed main_figure_beds/Sp9_2M.bed
cp PSWM_family_3_2M.bed main_figure_beds/Rorb_2M.bed
cp PSWM_family_4_2M.bed main_figure_beds/RUNX_2M.bed
cp PSWM_family_5_2M.bed main_figure_beds/TEAD4_2M.bed

cp PSWM_family_1_fimo_all.bed main_figure_beds/Bach2_fimo_all.bed
cp PSWM_family_2a_fimo_all.bed main_figure_beds/Klf1_fimo_all.bed
cp PSWM_family_2b_fimo_all.bed main_figure_beds/Sp9_fimo_all.bed
cp PSWM_family_3_fimo_all.bed main_figure_beds/Rorb_fimo_all.bed
cp PSWM_family_4_fimo_all.bed main_figure_beds/RUNX_fimo_all.bed
cp PSWM_family_5_fimo_all.bed main_figure_beds/TEAD4_fimo_all.bed
```

| PSWM_family | specific_factor | family_name |
|-------------|-----------------|-------------|
| 1           | Bach2           | bZip/AP-1   |
| 2a          | Klf1            | KLF         |
| 2b          | Sp9             | SP          |
| 3           | Rorb            | RORs        |
| 4           | RUNX            | RUNX        |
| 5           | TEAD4           | TEAD or TEF | 