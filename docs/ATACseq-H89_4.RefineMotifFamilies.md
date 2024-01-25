
### 4. Combine motifs of the same family manually

Several of the individual motif families should likely be combined since they are from the same family

| PSWM_family | specific_factor | family_name |
|-------------|-----------------|-------------|
| 1           | Fra2            | AP-1        |
| 2a          | Klf1            | KLF         |
| 2b          | Sp9             | SP          |
| 3           | Rorb            | RORs        |
| 4           | Esrrb           | NR3B        |
| 5           | ETS-RUNX        | ETS/RUNX    | -
| 6           | Nr5a2           | NR5A        |
| 7           | PU.1            | ETS         | -
| 8           | RXR             | NR2B        |
| 9           | Zfp281          | ZFPs        | * or are ZFPs, ZNFs (ZFPS are huge family). Also, KLFs are SFPs...
| 10          | ZNF467          | ZFPs        | *
| 11          | Patz1           | ZBTB        |
| 12          | Zfp422          | ZFPs        | *
| 13          | Zfp467          | ZFPs        | *

#### ETS

Combine family 5 and 7

```console
cd FIMO/tomtom/
cat PSWM_family_5.txt PSWM_family_7.txt > PSWM_ETS_family.txt
```

Generate composite motifs
 
```console
mkdir -p FIMO/tomtom/composite_motifs/PSWM_family_ETS/
cd FIMO/tomtom/composite_motifs/PSWM_family_ETS/
cat ../../PSWM_ETS_family.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_ETS_family_factors_meme.txt
```

Check first family member
```console
tomtom -no-ssc -oc 'ETS:RUNX_homer.tomtom_output' -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 'ETS:RUNX_homer_meme.txt' PSWM_ETS_family_factors_meme.txt && max_motif=$(wc -l < 'ETS:RUNX_homer.tomtom_output/tomtom.tsv') && final_query="ETS:RUNX_homer"
```

Check the rest of the family members in the `PSWM_ETS_family_factors_meme.txt` file
```console
tomtom -no-ssc -oc 'PU.1_homer.tomtom_output' -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 'PU.1_homer_meme.txt' PSWM_ETS_family_factors_meme.txt
if [[ $(wc -l < 'PU.1_homer.tomtom_output/tomtom.tsv') -ge $max_motif ]]
then
    max_motif=$(wc -l < 'PU.1_homer.tomtom_output/tomtom.tsv')
    final_query="PU.1_homer"
fi

echo FINAL_QUERY IS $final_query
```

FINAL_QUERY IS PU.1_homer (max_motif: 7)

```console
final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_ETS_family_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_ETS_family_composite_index.txt
cp ../PSWM_ETS_family_composite_values.txt ../composite.values.txt
cp ../PSWM_ETS_family_composite_index.txt ../composite.index.txt
``` 

#### ZFPs

family 9,10 and 12,13

```console
cd /FIMO/tomtom/
cat PSWM_family_9.txt PSWM_family_10.txt PSWM_family_12.txt PSWM_family_13.txt > PSWM_ZFPs_family.txt
```

Generate composite motifs
 
```console
mkdir -p FIMO/tomtom/composite_motifs/PSWM_family_ZFPs/
cd FIMO/tomtom/composite_motifs/PSWM_family_ZFPs/
cat ../../PSWM_ZFPs_family.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done
cat *_meme.txt > PSWM_ZFPs_family_factors_meme.txt
```

Check each family member, for example:

```console
tomtom -no-ssc -oc 'Zfp281_cisbp.tomtom_output' -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 'Zfp281_cisbp_meme.txt' PSWM_ZFPs_family_factors_meme.txt && max_motif=$(wc -l < 'Zfp281_cisbp.tomtom_output/tomtom.tsv') && final_query="Zfp281_cisbp"

wc -l < 'Zfp281_primary_uniprobe.tomtom_output/tomtom.tsv'
wc -l < 'Zfp740_primary_uniprobe.tomtom_output/tomtom.tsv'
wc -l < 'ZNF467_homer.tomtom_output/tomtom.tsv'
wc -l < '(Zfp182)_(Homo_sapiens)_(DBD_0.97)_cisbp.tomtom_output/tomtom.tsv'
wc -l < '(Zfp287)_(Homo_sapiens)_(DBD_0.93)_cisbp.tomtom_output/tomtom.tsv'
wc -l < '(Zfp384)_(Homo_sapiens)_(DBD_0.99)_cisbp.tomtom_output/tomtom.tsv'
wc -l < '(Zfp422)_(Homo_sapiens)_(DBD_0.93)_cisbp.tomtom_output/tomtom.tsv'
wc -l < '(Zfp467)_(Homo_sapiens)_(DBD_0.95)_cisbp.tomtom_output/tomtom.tsv'

echo FINAL_QUERY IS $final_query
```

FINAL_QUERY IS Zfp281_cisbp (max_motif: 9)

```console
final_query_dir="${final_query}.tomtom_output/"

wc -l ${final_query_dir}/tomtom.tsv
cd ${final_query_dir}
python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../PSWM_ZFPs_family_composite_values.txt
mv tomtom.xml_test_index_rc_offset.txt ../PSWM_ZFPs_family_composite_index.txt
cp ../PSWM_ZFPs_family_composite_values.txt ../composite.values.txt
cp ../PSWM_ZFPs_family_composite_index.txt ../composite.index.txt
``` 

 - Must update the `generate_composite_motif.R` script to have correct dir.

```console
cd FIMO/tomtom/composite_motifs/

for family in PSWM_family_[ZE][FT][PS]*
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

run FIMO against composite PSWM and take top 2 million hits
```console
cd FIMO/tomtom/fimo_composites
```

#### ZFPs Family

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_ZFPs/PSWM_family_ZFPs_meme.txt mm10.fa > PSWM_family_ZFPs_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_ZFPs_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_ZFPs_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_ZFPs_2M.txt

tomtom -no-ssc -oc PSWM_family_ZFPs_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_ZFPs_composite_fimo.txt
```

#### ETS Family

```console
fimo --thresh 0.001 --text ../composite_motifs/PSWM_family_ETS/PSWM_family_ETS_meme.txt mm10.fa > PSWM_family_ETS_composite_fimo.txt

#this takes top 2M
score=$(tail -n +2 PSWM_family_ETS_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
echo $score
tail -n +2 PSWM_family_ETS_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > PSWM_family_ETS_2M.txt

tomtom -no-ssc -oc PSWM_family_ETS_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt PSWM_family_ETS_composite_fimo.txt
```

```console
for i in PSWM_family_[ZE][FT][PS]_2M.txt
for i in PSWM_family_ZFPs_2M.txt
for i in PSWM_family_ETS_2M.txt
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_2M.txt" '{print $1}')
    echo $name
    intersectBed -loj -a ../../../dynamic_peaks_2023-07-03.bed -b $i > ${name}_fimo.bed
    intersectBed -loj -a ../../../nondynamic_peaks.bed -b $i > ${name}_fimo_nondyn.bed
    intersectBed -loj -a ../../../all_peaks.bed -b $i > ${name}_fimo_all.bed
    cat $i | cut -f1-3,5 | sort -k1,1 -k2,2n > ${name}_2M.bed 
done

cp PSWM_family_ZFPs_fimo.bed main_figure_beds/ZFPs_fimo.bed
cp PSWM_family_ETS_fimo.bed main_figure_beds/ETS_fimo.bed
cp PSWM_family_ZFPs_fimo_all.bed main_figure_beds/ZFPs_fimo_all.bed
cp PSWM_family_ETS_fimo_all.bed main_figure_beds/ETS_fimo_all.bed
cp PSWM_family_ZFPs_2M.bed main_figure_beds/ZFPs_2M.bed
cp PSWM_family_ETS_2M.bed main_figure_beds/ETS_2M.bed
```