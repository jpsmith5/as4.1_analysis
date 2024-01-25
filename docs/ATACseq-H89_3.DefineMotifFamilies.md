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
slopBed -i cluster_bed_cluster1.bed \
      -g mm10.chrom.sizes -b -50 | \
        fastaFromBed -fi mm10.fa -bed stdin \
              -fo cluster1.fasta
meme -p 16 -oc FIMO/cluster1_meme -nmotifs 15 -objfun classic -evt 0.01 -searchsize 0 -minw 6 \
         -maxw 18 -revcomp -dna -markov_order 3 -maxsize 100000000 \
         cluster1.fasta

## cluster2
slopBed -i cluster_bed_cluster2.bed \
      -g mm10.chrom.sizes -b -50 | \
        fastaFromBed -fi mm10.fa -bed stdin \
              -fo cluster2.fasta
meme -p 16 -oc FIMO/cluster2_meme -nmotifs 15 -objfun classic -evt 0.01 -searchsize 0 -minw 6 \
         -maxw 18 -revcomp -dna -markov_order 3 -maxsize 100000000 \
         cluster2.fasta
```

### 3b. Generate combined motif database

```console
mkdir databases/
cd databases/

#Generate homer_uniprobe_jaspar PSWM database
wget https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.23.tgz
tar -xzf motif_databases.12.23.tgz

#JASPAR
mv motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme $PWD
#Uniprobe
mv motif_databases/MOUSE/uniprobe_mouse.meme $PWD
#HOCOMOCO
mv motif_databases/MOUSE/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme $PWD
#CIS_BP
mv motif_databases/CIS-BP_2.00/Mus_musculus.meme $PWD

#Homer
#CAUTION: the HOMER_MEME_conversion.py was written for Python2 so remember to specify python2.7 when running.
wget https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/HOMER_MEME_conversion.py
wget http://homer.ucsd.edu/homer/custom.motifs
python2.7 HOMER_MEME_conversion.py -i custom.motifs -o homer.motifs

#edit databases to work with tomtom
## JASPAR
cp JASPAR2022_CORE_vertebrates_non-redundant_v2.meme JASPAR_edited_meme.txt
grep MOTIF JASPAR_edited_meme.txt > motifs.txt
cat motifs.txt | while read motif
do
    name=$(echo $motif | awk -F" " '{print $NF}')
    temp=$(echo 'MOTIF' $name'_jaspar')
    echo $temp
    sed -i "s;${motif};${temp};g" JASPAR_edited_meme.txt
done
rm motifs.txt

## uniprobe
cp uniprobe_mouse.meme uniprobe_edited_meme.txt
grep MOTIF uniprobe_edited_meme.txt > motifs.txt
cat motifs.txt | while read motif
do
    name=$(echo $motif | awk -F" " '{print $NF}')
    temp=$(echo 'MOTIF' $name'_uniprobe')
    echo $temp
    sed -i "s;${motif};${temp};g" uniprobe_edited_meme.txt
done
#remove 'secondary' motifs
grep -n 'secondary' uniprobe_edited_meme.txt
sed -i -e '4210,8326d;' uniprobe_edited_meme.txt

rm motifs.txt

## HOCOMOCO
cp HOCOMOCOv11_core_MOUSE_mono_meme_format.meme HOCOMOCO_edited_meme.txt
grep MOTIF HOCOMOCO_edited_meme.txt > motifs.txt
cat motifs.txt | while read motif
do
    name=$(echo $motif | awk -F" " '{print $NF}')
    temp=$(echo 'MOTIF' $name'_hocomoco')
    echo $temp
    sed -i "s;${motif};${temp};g" HOCOMOCO_edited_meme.txt
done
rm motifs.txt

## CIS_BP
cp Mus_musculus.meme CISBP_edited_meme.txt
grep MOTIF CISBP_edited_meme.txt > motifs.txt
cat motifs.txt | while read motif
do
    name=$(echo $motif | awk -F" " '{print $NF}')
    temp=$(echo 'MOTIF' $name'_cisbp')
    echo $temp
    sed -i "s;${motif};${temp};g" CISBP_edited_meme.txt
done

## HOMER
cp homer.motifs_meme.txt homer_edited_meme.txt
grep MOTIF homer_edited_meme.txt > motifs.txt
cat motifs.txt | while read motif
do
    name=$(echo $motif | awk -F" " '{print $NF}')
    name=$(echo $name | awk -F"/" '{print $1}')
    temp=$(echo 'MOTIF' $name'_homer')
    echo $temp
    sed -i "s;${motif};${temp};g" homer_edited_meme.txt
done
rm motifs.txt

#Collect all database motifs into one file
cat CISBP_edited_meme.txt homer_edited_meme.txt uniprobe_edited_meme.txt JASPAR_edited_meme.txt HOCOMOCO_edited_meme.txt > motif_database_combined.txt

cd ../
```

### 3c. Run TomTom on Meme output

```console
mkdir tomtom/
cd tomtom/

echo 'Running TOMTOM'

tomtom -no-ssc -oc cluster1.tomtom_output -verbosity 1 -min-overlap 5 \
 -dist ed -evalue -thresh 0.05 \
 ../cluster1_meme/meme.txt \
 ../databases/motif_database_combined.txt

tomtom -no-ssc -oc cluster2.tomtom_output -verbosity 1 -min-overlap 5 \
 -dist ed -evalue -thresh 0.05 \
 ../cluster2_meme/meme.txt \
 ../databases/motif_database_combined.txt

tomtom -no-ssc -oc cluster1.tomtom_output -verbosity 1 -min-overlap 5 \
 -dist ed -text -evalue -thresh 0.05 \
 ../cluster1_meme/meme.txt \
 ../databases/motif_database_combined.txt > cluster1.tomtom_output/tomtom.txt

tomtom -no-ssc -oc cluster2.tomtom_output -verbosity 1 -min-overlap 5 \
 -dist ed -text -evalue -thresh 0.05 \
 ../cluster2_meme/meme.txt \
 ../databases/motif_database_combined.txt > cluster2.tomtom_output/tomtom.txt
```

 - Extract motif names from tomtom output.
 - Extract individual meme files from combined database.

```console
mkdir individual_memes/
cd individual_memes/

wget https://github.com/guertinlab/adipogenesis/blob/2fdf0bbab4fe6f368f5a60e42f7899b6570ff71c/motif_analysis/MEME_individual_from_db.py
python2.7 MEME_individual_from_db.py -i ../../databases/motif_database_combined.txt

for file in *meme.txt 
do
    name=$(echo $file | awk -F"motif_database_combined.txt_" '{print $1}')
    mv $file ${name}meme.txt
    
done

#tomtom all query factors against all others
cd ../
mkdir tomtom_all_query_factors/
cd tomtom_all_query_factors/
```

`motif_names.txt` is a single column list of the motifs discovered from the above `tomtom` command. Make one for each identified cluster (i.e. cluster1 and cluster2).

Copy the ID column from tomtom.tsv: `cluster1_motif_names.txt`, `cluster2_motif_names.txt`

| cluster1                                 | cluster2                                 |
|------------------------------------------|------------------------------------------|
| AP-1_homer                               | Sp1_homer                                |
| Bach1_homer                              | Sp2_homer                                |
| Fra2_homer                               | (Klf15)_(Homo_sapiens)_(DBD_1.00)_cisbp  |
| Nrf2_homer                               | (Klf14)_(Homo_sapiens)_(DBD_0.97)_cisbp  |
| Atf3_homer                               | Sp5_homer                                |
| NF-E2_homer                              | (Sp2)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
| BATF_homer                               | KLF3_homer                               |
| (Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp | KLF1_homer                               |
| Jun-AP1_homer                            | KLF14_homer                              |
| (Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp    | (Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp  |
| Fos_homer                                | Klf7_primary_uniprobe                    |
| (Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp  | (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
| Bach2_homer                              | (Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp   |
| (Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp   | Zfp281_cisbp                             |
| Fra1_homer                               | Sp4_cisbp                                |
| JunB_homer                               | (Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp  |
| NFE2L2_homer                             | Sp4_primary_uniprobe                     |
| Jund_cisbp                               | Klf4_cisbp                               |
| Nfe2l2_cisbp                             | Zfp740_primary_uniprobe                  |
| Fosl2_homer                              | Sp1_cisbp                                |
| EKLF_homer                               | (Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp    |
| KLF5_homer                               | (Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp  |
| Klf1_cisbp                               | (Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
| Klf4_cisbp                               | (Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp |
| (Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp  | KLF5_homer                               |
| KLF6_homer                               | Klf7_cisbp                               |
| Klf4_homer                               | KLF6_homer                               |
| KLF1_homer                               | Zfp281_primary_uniprobe                  |
| (Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp   | (Klf11)_(Homo_sapiens)_(DBD_0.97)_cisbp  |
| KLF3_homer                               | (Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp   |
| (Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp   | Klf9_homer                               |
| (Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp    | Klf12_cisbp                              |
| (Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp  | (Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp   |
| Klf9_homer                               | ZNF467_homer                             |
| (Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp   | (Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp   |
| Klf7_primary_uniprobe                    | ERE_homer                                |
| (Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp    | Nr2c1_cisbp                              |
| Klf7_cisbp                               | Rxrb_cisbp                               |
| Klf12_cisbp                              | Rxrg_cisbp                               |
| (Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp   | Nr2f2_cisbp                              |
| KLF14_homer                              | Rorb_cisbp                               |
| KLF10_homer                              | Fra1_homer                               |
| (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp    | Fos_homer                                |
| ETS:RUNX_homer                           | Atf3_homer                               |
| PU.1_homer                               | Jun-AP1_homer                            |
|                                          | BATF_homer                               |
|                                          | Fra2_homer                               |
|                                          | (Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
|                                          | JunB_homer                               |
|                                          | Nfe2l2_cisbp                             |
|                                          | AP-1_homer                               |
|                                          | Fosl2_homer                              |
|                                          | (Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp |
|                                          | NF-E2_homer                              |
|                                          | (Batf)_(Homo_sapiens)_(DBD_0.98)_cisbp   |
|                                          | (Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp   |
|                                          | Bach2_homer                              |
|                                          | Bach1_homer                              |
|                                          | Nrf2_homer                               |
|                                          | (Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp  |
|                                          | NFAT:AP1_homer                           |
|                                          | NFE2L2_homer                             |
|                                          | (Zfp422)_(Homo_sapiens)_(DBD_0.93)_cisbp |
|                                          | (Zfp287)_(Homo_sapiens)_(DBD_0.93)_cisbp |
|                                          | (Zfp384)_(Homo_sapiens)_(DBD_0.99)_cisbp |
|                                          | (Zfp182)_(Homo_sapiens)_(DBD_0.97)_cisbp |
|                                          | Esrrb_homer                              |
|                                          | Nr5a2_homer                              |
|                                          | RXR_homer                                |

```console
cat ../cluster1_motif_names.txt | while read factor
do
    echo $factor
    cp ../individual_memes/${factor}_meme.txt $PWD
done

cat *meme.txt > ../cluster1_query_factors_meme.txt

cat ../cluster2_motif_names.txt | while read factor
do
    echo $factor
    cp ../individual_memes/${factor}_meme.txt $PWD
done

cat *meme.txt > ../cluster2_query_factors_meme.txt

cd ../
cat cluster1_query_factors_meme.txt cluster2_query_factors_meme.txt > all_query_factors_meme.txt
cd tomtom_all_query_factors/
```

### 3d. Define motif families

On any of the motifs with a parenthesis, may need to run `tomtom` manually. For example:

| cluster1                                 | cluster2                                 |
|------------------------------------------|------------------------------------------|
| (Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp  | (Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp  |
| (Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp    | (Batf)_(Homo_sapiens)_(DBD_0.98)_cisbp   |
| (Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp  | (Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
| (Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp  | (Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp  |
| (Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp   | (Klf11)_(Homo_sapiens)_(DBD_0.97)_cisbp  |
| (Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp   | (Klf14)_(Homo_sapiens)_(DBD_0.97)_cisbp  |
| (Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp   | (Klf15)_(Homo_sapiens)_(DBD_1.00)_cisbp  |
| (Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp   | (Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp  |
| (Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp   | (Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp   |
| (Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp | (Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp   |
| (Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp    | (Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp   |
| (Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp    | (Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp   |
| (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp    | (Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp   |
|                                          | (Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp |
|                                          | (Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp  |
|                                          | (Sp2)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
|                                          | (Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp    |
|                                          | (Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
|                                          | (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
|                                          | (Zfp182)_(Homo_sapiens)_(DBD_0.97)_cisbp |
|                                          | (Zfp287)_(Homo_sapiens)_(DBD_0.93)_cisbp |
|                                          | (Zfp384)_(Homo_sapiens)_(DBD_0.99)_cisbp |
|                                          | (Zfp422)_(Homo_sapiens)_(DBD_0.93)_cisbp |
|                                          | (Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp |

Otherwise, can run them in a loop.

```console
for meme in *meme.txt
do
    name=$(echo $meme | awk -F".txt_" '{print $NF}' | awk -F"_meme.txt" '{print $1}')
    echo $name
    tomtom -no-ssc -oc $name.tomtom_output -verbosity 1 \
     -incomplete-scores -min-overlap 1 -dist ed \
     -evalue -thresh 0.0005 ${meme} ../all_query_factors_meme.txt
    cd $name.tomtom_output
    cut -f1,2,5 tomtom.tsv | tail -n +2 | sed '$d' | sed '$d' \
     | sed '$d' | sed '$d' >> ../3_col_combined_motif_db_pre.txt
    cd ..
done

# May need to run tomtom on files with irregular characters. For example:
tomtom -no-ssc -oc '(Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp.tomtom_output' \
 -verbosity 1  -incomplete-scores -min-overlap 1 -dist ed -evalue \
 -thresh 0.0005 '(Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp_meme.txt' \
 ../all_query_factors_meme.txt
cd '(Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp.tomtom_output' && \
 cut -f1,2,5 tomtom.tsv | tail -n +2 | sed '$d' | sed '$d' | \
 sed '$d' | sed '$d' >> ../3_col_combined_motif_db_pre.txt && cd ../

grep -v '#' 3_col_combined_motif_db_pre.txt > 3_col_combined_motif_db.txt
rm 3_col_combined_motif_db_pre.txt
cp 3_col_combined_motif_db.txt ..
cd ..
```

Split the SP/KLF family into individual families.

```console
cd /scratch/jps3dp/processed/pepatac/gomez/robert/H89/bam_files/FIMO/tomtom/
grep -e '[Kk][Ll][Ff]' PSWM_family_2.txt > PSWM_family_2a.txt
grep -e '[Ss][Pp]' PSWM_family_2.txt > PSWM_family_2b.txt
```

| Family1                                  | Family2                                 | Family2a                                | Family2b                              | Family3     | Family4     | Family5        | Family6     | Family7    | Family8   | Family9                 | Family10     | Family11                                | Family12                                 | Family13                                 |
|------------------------------------------|-----------------------------------------|-----------------------------------------|---------------------------------------|-------------|-------------|----------------|-------------|------------|-----------|-------------------------|--------------|-----------------------------------------|------------------------------------------|------------------------------------------|
| AP-1_homer                               | EKLF_homer                              | EKLF_homer                              | Sp1_cisbp                             | ERE_homer   | Esrrb_homer | ETS:RUNX_homer | Nr5a2_homer | PU.1_homer | RXR_homer | Zfp281_cisbp            | ZNF467_homer | (Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp | (Zfp182)_(Homo_sapiens)_(DBD_0.97)_cisbp | (Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp |
| Atf3_homer                               | KLF10_homer                             | KLF10_homer                             | Sp1_homer                             | Nr2c1_cisbp |             |                |             |            |           | Zfp281_primary_uniprobe |              |                                         | (Zfp287)_(Homo_sapiens)_(DBD_0.93)_cisbp |                                          |
| Bach1_homer                              | Klf12_cisbp                             | Klf12_cisbp                             | Sp2_homer                             | Nr2f2_cisbp |             |                |             |            |           | Zfp740_primary_uniprobe |              |                                         | (Zfp384)_(Homo_sapiens)_(DBD_0.99)_cisbp |                                          |
| Bach2_homer                              | KLF14_homer                             | KLF14_homer                             | Sp4_cisbp                             | Rorb_cisbp  |             |                |             |            |           |                         |              |                                         | (Zfp422)_(Homo_sapiens)_(DBD_0.93)_cisbp |                                          |
| BATF_homer                               | Klf1_cisbp                              | Klf1_cisbp                              | Sp4_primary_uniprobe                  | Rxrb_cisbp  |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| Fos_homer                                | KLF1_homer                              | KLF1_homer                              | Sp5_homer                             | Rxrg_cisbp  |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| Fosl2_homer                              | KLF3_homer                              | KLF3_homer                              | (Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| Fra1_homer                               | Klf4_cisbp                              | Klf4_cisbp                              | (Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| Fra2_homer                               | Klf4_homer                              | Klf4_homer                              | (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| Jun-AP1_homer                            | KLF5_homer                              | KLF5_homer                              | (Sp2)_(Homo_sapiens)_(DBD_1.00)_cisbp |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| JunB_homer                               | KLF6_homer                              | KLF6_homer                              |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| Jund_cisbp                               | Klf7_cisbp                              | Klf7_cisbp                              |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| NFAT:AP1_homer                           | Klf7_primary_uniprobe                   | Klf7_primary_uniprobe                   |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| NF-E2_homer                              | Klf9_homer                              | Klf9_homer                              |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| Nfe2l2_cisbp                             | Sp1_cisbp                               | (Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| NFE2L2_homer                             | Sp1_homer                               | (Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| Nrf2_homer                               | Sp2_homer                               | (Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp  |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| (Bach2)_(Homo_sapiens)_(DBD_0.95)_cisbp  | Sp4_cisbp                               | (Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp  |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| (Fos)_(Homo_sapiens)_(DBD_1.00)_cisbp    | Sp4_primary_uniprobe                    | (Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp  |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| (Nfe2)_(Homo_sapiens)_(DBD_0.95)_cisbp   | Sp5_homer                               | (Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp  |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| (Nfe2l1)_(Homo_sapiens)_(DBD_1.00)_cisbp | (Klf10)_(Homo_sapiens)_(DBD_0.99)_cisbp | (Klf11)_(Homo_sapiens)_(DBD_0.97)_cisbp |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
| (Batf)_(Homo_sapiens)_(DBD_0.98)_cisbp   | (Klf16)_(Homo_sapiens)_(DBD_0.94)_cisbp | (Klf14)_(Homo_sapiens)_(DBD_0.97)_cisbp |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Klf2)_(Homo_sapiens)_(DBD_0.99)_cisbp  | (Klf15)_(Homo_sapiens)_(DBD_1.00)_cisbp |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Klf3)_(Homo_sapiens)_(DBD_1.00)_cisbp  |                                         |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Klf5)_(Homo_sapiens)_(DBD_0.96)_cisbp  |                                         |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Klf6)_(Homo_sapiens)_(DBD_1.00)_cisbp  |                                         |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Sp3)_(Homo_sapiens)_(DBD_0.97)_cisbp   |                                         |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Sp8)_(Homo_sapiens)_(DBD_1.00)_cisbp   |                                         |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp   |                                         |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Klf11)_(Homo_sapiens)_(DBD_0.97)_cisbp |                                         |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Klf14)_(Homo_sapiens)_(DBD_0.97)_cisbp |                                         |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Klf15)_(Homo_sapiens)_(DBD_1.00)_cisbp |                                         |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |
|                                          | (Sp2)_(Homo_sapiens)_(DBD_1.00)_cisbp   |                                         |                                       |             |             |                |             |            |           |                         |              |                                         |                                          |                                          |

```console
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Pipeline_ATAC/misc_scripts/tomtom_output_to_composite.py
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Pipeline_ATAC/misc_scripts/meme_header.txt
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Pipeline_ATAC/misc_scripts/generate_composite_motif.R

mkdir composite_motifs/
cd composite_motifs/

mkdir PSWM_family_1/
mkdir PSWM_family_2/
mkdir PSWM_family_2a/
mkdir PSWM_family_2b/
mkdir PSWM_family_3/
mkdir PSWM_family_4/
mkdir PSWM_family_5/
mkdir PSWM_family_6/
mkdir PSWM_family_7/
mkdir PSWM_family_8/
mkdir PSWM_family_9/
mkdir PSWM_family_10/
mkdir PSWM_family_11/
mkdir PSWM_family_12/
mkdir PSWM_family_13/
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

# Check first family member (then repeat for each additional member)
tomtom -no-ssc -oc AP-1_homer.tomtom_output -verbosity 1 -incomplete-scores \
 -min-overlap 1 -dist ed -evalue -thresh 0.0005 \
 AP-1_homer_meme.txt PSWM_family_1_factors_meme.txt

max_motif=$(wc -l < AP-1_homer.tomtom_output/tomtom.tsv)
final_query="AP-1_homer"

# Check the rest ...
echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS: Fra2_homer

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_1_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_1_composite_index.txt
```

Now, for each additional family, repeat the same process.

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
tomtom -no-ssc -oc 'EKLF_homer.tomtom_output' -verbosity 1 -incomplete-scores \
 -min-overlap 1 -dist ed -evalue -thresh 0.0005 'EKLF_homer_meme.txt' \
 PSWM_family_2_factors_meme.txt && max_motif=$(wc -l < \
 'EKLF_homer.tomtom_output/tomtom.tsv') && final_query="EKLF_homer"

# Check the rest ...

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS Klf1_cisbp

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_2_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_2_composite_index.txt
``` 

##### KLF (2a)
```console
cd ../../PSWM_family_2a/
cat ../../PSWM_family_2a.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_2a_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'EKLF_homer.tomtom_output' -verbosity 1 -incomplete-scores \
 -min-overlap 1 -dist ed -evalue -thresh 0.0005 'EKLF_homer_meme.txt' \
  PSWM_family_2a_factors_meme.txt && max_motif=$(wc -l < \
  'EKLF_homer.tomtom_output/tomtom.tsv') && final_query="EKLF_homer"

# Check the rest ...

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS Klf1_cisbp (max_motif: 20)

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_2a_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_2a_composite_index.txt
cp ../PSWM_family_2a_composite_values.txt ../composite.values.txt
cp ../PSWM_family_2a_composite_index.txt ../composite.index.txt
``` 

##### SP (2b)
```console
cd ../../PSWM_family_2b/
cat ../../PSWM_family_2b.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_2b_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'Sp1_cisbp.tomtom_output' -verbosity 1 -incomplete-scores \
 -min-overlap 1 -dist ed -evalue -thresh 0.0005 'Sp1_cisbp_meme.txt' \
  PSWM_family_2b_factors_meme.txt && max_motif=$(wc -l < \
  'Sp1_cisbp.tomtom_output/tomtom.tsv') && final_query="Sp1_cisbp"

# Check the rest ...

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp (max_motif: 9)

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_2b_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_2b_composite_index.txt
cp ../PSWM_family_2b_composite_values.txt ../composite.values.txt
cp ../PSWM_family_2b_composite_index.txt ../composite.index.txt
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
tomtom -no-ssc -oc 'ERE_homer.tomtom_output' -verbosity 1 -incomplete-scores \
 -min-overlap 1 -dist ed -evalue -thresh 0.0005 'ERE_homer_meme.txt' \
 PSWM_family_3_factors_meme.txt && max_motif=$(wc -l < \
 'ERE_homer.tomtom_output/tomtom.tsv') && final_query="ERE_homer"

# Check the rest ...

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS Rorb_cisbp

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
tomtom -no-ssc -oc 'Esrrb_homer.tomtom_output' -verbosity 1 -incomplete-scores \
 -min-overlap 1 -dist ed -evalue -thresh 0.0005 'Esrrb_homer_meme.txt' \
 PSWM_family_4_factors_meme.txt && max_motif=$(wc -l < \
 'Esrrb_homer.tomtom_output/tomtom.tsv') && final_query="Esrrb_homer"

# There's only 1

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS Esrrb_homer

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
tomtom -no-ssc -oc 'ETS:RUNX_homer.tomtom_output' -verbosity 1 \
 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 \
 'ETS:RUNX_homer_meme.txt' PSWM_family_5_factors_meme.txt && \
 max_motif=$(wc -l < 'ETS:RUNX_homer.tomtom_output/tomtom.tsv') && \
 final_query="ETS:RUNX_homer"


# There's only 1

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS ETS:RUNX_homer

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_5_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_5_composite_index.txt
``` 

#### PSWM Family 6
```console
cd ../../PSWM_family_6/
cat ../../PSWM_family_6.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_6_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'Nr5a2_homer.tomtom_output' -verbosity 1 -incomplete-scores \
 -min-overlap 1 -dist ed -evalue -thresh 0.0005 'Nr5a2_homer_meme.txt' \
 PSWM_family_6_factors_meme.txt && max_motif=$(wc -l < \
 'Nr5a2_homer.tomtom_output/tomtom.tsv') && final_query="Nr5a2_homer"

# There's only 1

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS Nr5a2_homer

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_6_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_6_composite_index.txt
``` 

#### PSWM Family 7
```console
cd ../../PSWM_family_7/
cat ../../PSWM_family_7.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_7_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'PU.1_homer.tomtom_output' -verbosity 1 -incomplete-scores \
 -min-overlap 1 -dist ed -evalue -thresh 0.0005 'PU.1_homer_meme.txt' \
 PSWM_family_7_factors_meme.txt && max_motif=$(wc -l < \
 'PU.1_homer.tomtom_output/tomtom.tsv') && final_query="PU.1_homer"

# There's only 1

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS PU.1_homer

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_7_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_7_composite_index.txt
``` 

#### PSWM Family 8
```console
cd ../../PSWM_family_8/
cat ../../PSWM_family_8.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_8_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'RXR_homer.tomtom_output' -verbosity 1 -incomplete-scores \
 -min-overlap 1 -dist ed -evalue -thresh 0.0005 'RXR_homer_meme.txt' \
 PSWM_family_8_factors_meme.txt && max_motif=$(wc -l < \
 'RXR_homer.tomtom_output/tomtom.tsv') && final_query="RXR_homer"

# There's only 1

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS RXR_homer

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_8_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_8_composite_index.txt
``` 

#### PSWM Family 9
```console
cd ../../PSWM_family_9/
cat ../../PSWM_family_9.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_9_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'Zfp281_cisbp.tomtom_output' -verbosity 1 \
 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 \
 'Zfp281_cisbp_meme.txt' PSWM_family_9_factors_meme.txt && \
 max_motif=$(wc -l < 'Zfp281_cisbp.tomtom_output/tomtom.tsv') && \
 final_query="Zfp281_cisbp"

# Check the rest ...

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS Zfp281_primary_uniprobe

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_9_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_9_composite_index.txt
```

#### PSWM Family 10
```console
cd ../../PSWM_family_10/
cat ../../PSWM_family_10.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_10_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc 'ZNF467_homer.tomtom_output' -verbosity 1 \
 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 \
 'ZNF467_homer_meme.txt' PSWM_family_10_factors_meme.txt && \
 max_motif=$(wc -l < 'ZNF467_homer.tomtom_output/tomtom.tsv') && \
 final_query="ZNF467_homer"

# Only 1

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS ZNF467_homer

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_10_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_10_composite_index.txt
```

#### PSWM Family 11
```console
cd ../../PSWM_family_11/
cat ../../PSWM_family_11.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_11_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc '(Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp.tomtom_output' \
 -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue \
 -thresh 0.0005 '(Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp_meme.txt' \
 PSWM_family_11_factors_meme.txt && max_motif=$(wc -l < \
 '(Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp.tomtom_output/tomtom.tsv') && \
 final_query="(Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp"

# Only 1

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS (Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_11_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_11_composite_index.txt
``` 

#### PSWM Family 12
```console
cd ../../PSWM_family_12/
cat ../../PSWM_family_12.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_12_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc '(Zfp182)_(Homo_sapiens)_(DBD_0.97)_cisbp.tomtom_output' \
 -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue \
 -thresh 0.0005 '(Zfp182)_(Homo_sapiens)_(DBD_0.97)_cisbp_meme.txt' \
 PSWM_family_12_factors_meme.txt && max_motif=$(wc -l < \
 '(Zfp182)_(Homo_sapiens)_(DBD_0.97)_cisbp.tomtom_output/tomtom.tsv') && \
 final_query="(Zfp182)_(Homo_sapiens)_(DBD_0.97)_cisbp"

# Only 1

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS (Zfp422)_(Homo_sapiens)_(DBD_0.93)_cisbp

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_12_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_12_composite_index.txt
```

#### PSWM Family 13
```console
cd ../../PSWM_family_13/
cat ../../PSWM_family_13.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_family_13_factors_meme.txt

# Check first family member
tomtom -no-ssc -oc '(Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp.tomtom_output' \
 -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue \
 -thresh 0.0005 '(Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp_meme.txt' \
 PSWM_family_13_factors_meme.txt && max_motif=$(wc -l < \
 '(Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp.tomtom_output/tomtom.tsv') && \
 final_query="(Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp"

# Only 1

echo FINAL_QUERY IS $final_query
# FINAL_QUERY IS (Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp

final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_family_13_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_family_13_composite_index.txt
```

```console
cd FIMO/tomtom/composite_motifs/

cp PSWM_family_1/PSWM_family_1_composite_index.txt PSWM_family_1/composite.index.txt
cp PSWM_family_1/PSWM_family_1_composite_values.txt PSWM_family_1/composite.values.txt        
cp PSWM_family_2/PSWM_family_2_composite_index.txt PSWM_family_2/composite.index.txt
cp PSWM_family_2/PSWM_family_2_composite_values.txt PSWM_family_2/composite.values.txt
cp PSWM_family_3/PSWM_family_3_composite_values.txt PSWM_family_3/composite.values.txt
cp PSWM_family_3/PSWM_family_3_composite_index.txt PSWM_family_3/composite.index.txt
cp PSWM_family_4/PSWM_family_4_composite_values.txt PSWM_family_4/composite.values.txt
cp PSWM_family_4/PSWM_family_4_composite_index.txt PSWM_family_4/composite.index.txt
cp PSWM_family_5/PSWM_family_5_composite_values.txt PSWM_family_5/composite.values.txt
cp PSWM_family_5/PSWM_family_5_composite_index.txt PSWM_family_5/composite.index.txt
cp PSWM_family_6/PSWM_family_6_composite_values.txt PSWM_family_6/composite.values.txt
cp PSWM_family_6/PSWM_family_6_composite_index.txt PSWM_family_6/composite.index.txt
cp PSWM_family_7/PSWM_family_7_composite_values.txt PSWM_family_7/composite.values.txt
cp PSWM_family_7/PSWM_family_7_composite_index.txt PSWM_family_7/composite.index.txt
cp PSWM_family_8/PSWM_family_8_composite_values.txt PSWM_family_8/composite.values.txt
cp PSWM_family_8/PSWM_family_8_composite_index.txt PSWM_family_8/composite.index.txt
cp PSWM_family_9/PSWM_family_9_composite_values.txt PSWM_family_9/composite.values.txt
cp PSWM_family_9/PSWM_family_9_composite_index.txt PSWM_family_9/composite.index.txt
cp PSWM_family_10/PSWM_family_10_composite_values.txt PSWM_family_10/composite.values.txt
cp PSWM_family_10/PSWM_family_10_composite_index.txt PSWM_family_10/composite.index.txt
cp PSWM_family_11/PSWM_family_11_composite_values.txt PSWM_family_11/composite.values.txt
cp PSWM_family_11/PSWM_family_11_composite_index.txt PSWM_family_11/composite.index.txt
cp PSWM_family_12/PSWM_family_12_composite_values.txt PSWM_family_12/composite.values.txt
cp PSWM_family_12/PSWM_family_12_composite_index.txt PSWM_family_12/composite.index.txt
cp PSWM_family_13/PSWM_family_13_composite_values.txt PSWM_family_13/composite.values.txt
cp PSWM_family_13/PSWM_family_13_composite_index.txt PSWM_family_13/composite.index.txt

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
| 1           | Fra2_homer                               |
| 2           | Klf1_cisbp                               |
| 2a          | Klf1_cisbp                               |
| 2b          | (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
| 3           | Rorb_cisbp                               |
| 4           | Esrrb_homer                              |
| 5           | ETS:RUNX_homer                           |
| 6           | Nr5a2_homer                              |
| 7           | PU.1_homer                               |
| 8           | RXR_homer                                |
| 9           | Zfp281_primary_uniprobe                  |
| 10          | ZNF467_homer                             |
| 11          | (Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp  |
| 12          | (Zfp422)_(Homo_sapiens)_(DBD_0.93)_cisbp |
| 13          | (Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp |

### 3e. Run FIMO against composite PSWM and take top 2 million hits

```console
cd FIMO/tomtom
mkdir fimo_composites/
cd fimo_composites/
```

Family 1 (Fra2_homer)

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_1/PSWM_family_1_meme.txt mm10.fa > PSWM_family_1_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_1_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_1_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_1_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_1_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_1_composite_fimo.txt
```

Family 2 (Klf1_cisbp)

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_2/PSWM_family_2_meme.txt mm10.fa > PSWM_family_2_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_2_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_2_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_2_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_2_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_2_composite_fimo.txt
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

Family 3

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_3/PSWM_family_3_meme.txt mm10.fa > PSWM_family_3_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_3_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_3_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_3_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_3_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_3_composite_fimo.txt
```

Family 4

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_4/PSWM_family_4_meme.txt mm10.fa > PSWM_family_4_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_4_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_4_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_4_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_4_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_4_composite_fimo.txt
```

Family 5

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_5/PSWM_family_5_meme.txt mm10.fa > PSWM_family_5_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_5_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_5_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_5_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_5_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_5_composite_fimo.txt
```

Family 6

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_6/PSWM_family_6_meme.txt mm10.fa > PSWM_family_6_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_6_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_6_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_6_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_6_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_6_composite_fimo.txt
```

Family 7

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_7/PSWM_family_7_meme.txt mm10.fa > PSWM_family_7_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_7_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_7_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_7_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_7_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_7_composite_fimo.txt
```

Family 8

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_8/PSWM_family_8_meme.txt mm10.fa > PSWM_family_8_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_8_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_8_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_8_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_8_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_8_composite_fimo.txt

# The above is claiming undefined alphabet, not clear why. Create BED manually
awk -F'\t' '{print $1, $2, $3, $5}' PSWM_family_8_2M.txt > PSWM_family_8_2M.bed
```

Family 9

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_9/PSWM_family_9_meme.txt mm10.fa > PSWM_family_9_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_9_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_9_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_9_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_9_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_9_composite_fimo.txt
```

Family 10

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_10/PSWM_family_10_meme.txt mm10.fa > PSWM_family_10_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_10_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_10_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_10_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_10_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_10_composite_fimo.txt
```

Family 11

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_11/PSWM_family_11_meme.txt mm10.fa > PSWM_family_11_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_11_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_11_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_11_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_11_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_11_composite_fimo.txt
```

Family 12

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_12/PSWM_family_12_meme.txt mm10.fa > PSWM_family_12_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_12_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_12_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_12_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_12_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_12_composite_fimo.txt
```

Family 13

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_13/PSWM_family_13_meme.txt mm10.fa > PSWM_family_13_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_13_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_13_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_13_2M.txt

#this was to get the order of conformity to consensus.
tomtom -no-ssc -oc PSWM_family_13_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_13_composite_fimo.txt
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
| 1           | Fra2_homer                               |
| 2           | Klf1_cisbp                               |
| 2a          | Klf1_cisbp                               |
| 2b          | (Sp9)_(Homo_sapiens)_(DBD_1.00)_cisbp    |
| 3           | Rorb_cisbp                               |
| 4           | Esrrb_homer                              |
| 5           | ETS:RUNX_homer                           |
| 6           | Nr5a2_homer                              |
| 7           | PU.1_homer                               |
| 8           | RXR_homer                                |
| 9           | Zfp281_primary_uniprobe                  |
| 10          | ZNF467_homer                             |
| 11          | (Patz1)_(Homo_sapiens)_(DBD_1.00)_cisbp  |
| 12          | (Zfp422)_(Homo_sapiens)_(DBD_0.93)_cisbp |
| 13          | (Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp |

Copy BED files into figure folder

```console
cp PSWM_family_1_fimo.bed main_figure_beds/Fra2_fimo.bed
cp PSWM_family_2a_fimo.bed main_figure_beds/Klf1_fimo.bed
cp PSWM_family_2b_fimo.bed main_figure_beds/Sp9_fimo.bed
cp PSWM_family_3_fimo.bed main_figure_beds/Rorb_fimo.bed
cp PSWM_family_4_fimo.bed main_figure_beds/Esrrb_fimo.bed
cp PSWM_family_5_fimo.bed main_figure_beds/ETS-RUNX_fimo.bed
cp PSWM_family_6_fimo.bed main_figure_beds/Nr5a2_fimo.bed
cp PSWM_family_7_fimo.bed main_figure_beds/PU.1_fimo.bed
cp PSWM_family_8_fimo.bed main_figure_beds/RXR_fimo.bed
cp PSWM_family_9_fimo.bed main_figure_beds/Zfp281_fimo.bed
cp PSWM_family_10_fimo.bed main_figure_beds/ZNF467_fimo.bed
cp PSWM_family_11_fimo.bed main_figure_beds/Patz1_fimo.bed
cp PSWM_family_12_fimo.bed main_figure_beds/Zfp422_fimo.bed
cp PSWM_family_13_fimo.bed main_figure_beds/Zfp467_fimo.bed

cp PSWM_family_1_2M.bed main_figure_beds/Fra2_2M.bed
cp PSWM_family_2a_2M.bed main_figure_beds/Klf1_2M.bed
cp PSWM_family_2b_2M.bed main_figure_beds/Sp9_2M.bed
cp PSWM_family_3_2M.bed main_figure_beds/Rorb_2M.bed
cp PSWM_family_4_2M.bed main_figure_beds/Esrrb_2M.bed
cp PSWM_family_5_2M.bed main_figure_beds/ETS-RUNX_2M.bed
cp PSWM_family_6_2M.bed main_figure_beds/Nr5a2_2M.bed
cp PSWM_family_7_2M.bed main_figure_beds/PU.1_2M.bed
cp PSWM_family_8_2M.bed main_figure_beds/RXR_2M.bed
cp PSWM_family_9_2M.bed main_figure_beds/Zfp281_2M.bed
cp PSWM_family_10_2M.bed main_figure_beds/ZNF467_2M.bed
cp PSWM_family_11_2M.bed main_figure_beds/Patz1_2M.bed
cp PSWM_family_12_2M.bed main_figure_beds/Zfp422_2M.bed
cp PSWM_family_13_2M.bed main_figure_beds/Zfp467_2M.bed

#rename 2M files for generating final ATAC dataframe
cp PSWM_family_1_fimo_all.bed main_figure_beds/Fra2_fimo_all.bed
cp PSWM_family_2a_fimo_all.bed main_figure_beds/Klf1_fimo_all.bed
cp PSWM_family_2b_fimo_all.bed main_figure_beds/Sp9_fimo_all.bed
cp PSWM_family_3_fimo_all.bed main_figure_beds/Rorb_fimo_all.bed
cp PSWM_family_4_fimo_all.bed main_figure_beds/Esrrb_fimo_all.bed
cp PSWM_family_5_fimo_all.bed main_figure_beds/ETS-RUNX_fimo_all.bed
cp PSWM_family_6_fimo_all.bed main_figure_beds/Nr5a2_fimo_all.bed
cp PSWM_family_7_fimo_all.bed main_figure_beds/PU.1_fimo_all.bed
cp PSWM_family_8_fimo_all.bed main_figure_beds/RXR_fimo_all.bed
cp PSWM_family_9_fimo_all.bed main_figure_beds/Zfp281_fimo_all.bed
cp PSWM_family_10_fimo_all.bed main_figure_beds/ZNF467_fimo_all.bed
cp PSWM_family_11_fimo_all.bed main_figure_beds/Patz1_fimo_all.bed
cp PSWM_family_12_fimo_all.bed main_figure_beds/Zfp422_fimo_all.bed
cp PSWM_family_13_fimo_all.bed main_figure_beds/Zfp467_fimo_all.bed
```