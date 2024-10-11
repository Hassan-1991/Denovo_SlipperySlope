#CLUSTERING METHODS#

#The following method clusters 1,907,201 proteins
#Of these, 634 are assigned to 2 clusters, while the rest just one cluster
#This leaves 141,621 proteins that are found in the 419 gff files, unassigned
#This could be because another protein from the same genome has been included in the cluster
#Or because searches weren't made against genomes which harbor the gene
#Using a 80/60/0.001 criteria, a further 27,065 proteins are assigned to unique clusters
grep -f proteins_missing_in_clusters.tsv 30547_vs_419proteome_ultrasensitive.tsv | awk -F '\t' '($5>60&&$3>80&&$16<0.001)' | cut -f1,2 | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep -v "," | awk -F '\t' '{OFS=FS}{print $2,$1}' > 27065_secondstep_clustering.tsv
#While 74,554 are assigned to >1 cluster
grep -f proteins_missing_in_clusters.tsv 30547_vs_419proteome_ultrasensitive.tsv | awk -F '\t' '($5>60&&$3>80&&$16<0.001)' | cut -f1,2 | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep "," | cut -f1 | sort -u | grep -f - 30547_vs_419proteome_ultrasensitive.tsv | awk -F '\t' '($5>60&&$3>80&&$16<0.001)' | cut -f1,2 > 74554_secondstep_clustering.tsv
#This leaves 40,002 unassigned, for which mmseqs/linclust was run (80/60 mode), finding 7214 unique sequences
mmseqs createdb 40002_nocluster.faa 40002_nocluster
mmseqs search 40002_nocluster 40002_nocluster resultDB tmp --min-seq-id 0.8 -c 0.6 --cov-mode 1
mmseqs convertalis 40002_nocluster 40002_nocluster resultDB resultDB.m8
mmseqs linclust 40002_nocluster clusterDB tmp --min-seq-id 0.8 -c 0.6 --cov-mode 1
mmseqs createtsv 40002_nocluster 40002_nocluster clusterDB 40002_clusters.tsv
awk '{if ($1 != prev) {count++; prev=$1}; print count, $0}' 40002_clusters.tsv | sed "s/^/mmseqs_/g" | sed "s/ /\t/g" > 40002_clusters.reps.tsv
#Final list of gene-cluster assignments:
cut -f1,3 40002_clusters.reps.tsv | cat - 74554_secondstep_clustering.tsv 27065_secondstep_clustering.tsv protein_clusters.tsv | cut -f-2 > 041124_all_protein_clusters.tsv
#For the ORFan analysis, all I need to is search the 7214 new clusters/representative sequences
#Methodology for generating a (i) pangenome, and a (ii) presence/absence table for E. coli:

"""1. Download 500 genomes used in the ipoppunk paper
-   I require lineage information that's why
2. Of these 500, 31 were missing from the presence/absence table that came with the OG paper
-   I wanna rely on their presence/absence data that's why
3. Of the remaining 461, 42 were the only representative genome in their respective lineage
-   I want to identify core genes in lineages vs. out, so let's get rid of them
4. Of the remaining 419, I identified genes in each cluster of representative genomes in the following way:
a) First, do a 90/90 search with the rep. genome as query and all prodigal-identified genes in the "present" genomes as target
-   In case of paralogs, pick the best hit
b) Second, do a 80/80 search with the rep. genome as query and all prodigal-identified genes in the "present" genomes as target
-   In case of paralogs, pick the best hit
c) Finally, do an 80/60 search
-   Compare how many strains are missing: if not more than 20%, then these diamond results are actual presence/absence results for this last third of genes"""

#Starting with the given presence/absence table
#Only retain those from the 500 who are present in the p/a table
#Only retain those which have 9/10 genomes from the same lineage

head -1 F4_complete_presence_absence.csv | sed "s/,/\n/g" | tail -n+2 | sed "s/$/\t/g" | grep -f - 500_ipp_lineagedesignations.tsv | awk -F '\t' '($5<45)' > 419_ipp_lineagedesignations.tsv
cut -f1 419_ipp_lineagedesignations.tsv | sed "1s/^/Strain\n/g" > 419_isolate_names.txt

#Take their p/a data
python3 parsing1.py 419_isolate_names.txt F4_complete_presence_absence.csv > 419_isolates_presence_absence.csv
#Remove all cases where the gene is absent from the entire 419-genome pangenome
grep ",1" 419_isolates_presence_absence.csv > 419_isolates_presence_absence.blanksremoved.csv
#Modify for further processing
sed -i "1s/,/,Strain_/g" 419_isolates_presence_absence.blanksremoved.csv
python3 parsing2.py 419_isolates_presence_absence.blanksremoved.csv > 419_isolates_presence_absence.blanksremoved.mod.csv
#parsing1.py and parsing2.py - paste the scripts here later - in an appendix?

#Overall strategy: search the genomes where the rep sequence is supposed to be present to identify its homolog

#First, make a list of files that list the genomes where each protein is supposed to be found:

tail -n+3 419_isolates_presence_absence.blanksremoved.mod.csv | head -n-1 |
awk -F',' '{for(i=2; i<=NF; i++) print $1","$i}' | #Convert file to long format
awk -F ',' '($2!=0)' | #Only retain presences - information about which isolates the gene is found in
sed "s/,/\t/g" | #convert to tab-separated
sed "s/\tStrain_/\t/g" > all_clusters_per_pa_table.tsv #Remove the "Strain" string we used earlier
#Run the following to see if you get the same number of lines as this tsv file
tail -n+3 419_isolates_presence_absence.blanksremoved.mod.csv | head -n-1 | cut -f2- -d ',' | sed "s/,/\n/g" | grep -v "^0" | wc -l

#Move this file to another directory
cp all_clusters_per_pa_table.tsv validation_positive/
cd validation_positive/

#Convert this master file into gene-by-gene files, where each file documents the "target" isolates for a particular gene
#One gene name has a "/", this needs to be modified everywhere protein names are listed
sed -i "s/rbsK\/rbiA/rbsK%rbiA/g" all_clusters_per_pa_table.tsv
awk -F'\t' '{print > $1".tsv"}' all_clusters_per_pa_table.tsv #Generates the per-gene files
for i in $(ls *tsv | grep -v "cluster"); do sort -nk 2 $i -o $i; done
cd ..

#This generates 30,547 files, one for each gene in the EC419 pangenome
#This directory can remain undisturbed

#Now for the searches:
#Search all representative protein sequences against the genomes extracted from 500 gff files
#This returns the file 500_gffs/500_proteins/ref_vs_500concat.tsv
#In cases where it has >1 hit, we retain the best hit (highest bit score)

#Strategy:
a) 80/60, only pick cases without paralogs
b) In case of paralogs, only pick the best hit
In case of a and b, we're only picking instances where the diamond-distribution exactly matches the p/a table distribution
c) For the remainder, pick the best hit, and if the number of strains detected is fewer than 10% of the total, it should be fine - pick the diamond p/a in these cases

for i in $(cut -f1 -d ',' 419_isolates_presence_absence.blanksremoved.mod.csv | tail -n+3 | head -n-1 | sed "s/$/,/g" | sed "s/,$//g")
do
echo "awk -v var=\"$i\" -F '\t' '(\$1==var&&\$3>80&&\$5>60)' 500_gffs/500_proteins/ref_vs_500concat.tsv | sed \"s/\t/@/\" | awk -F'@' 'NR==FNR {targets[\$0]; next} \$2 in targets' targets_in_pa_table/${i},@targets - | cut -f1 | cut -f2- -d \"@\" | sed \"s/^/${i}\t/g\" >> clusters.tsv"; done | 
awk '{print $0, int((NR-1)/25)+1}' | sed "s/clusters.tsv /clusters.tsv_/g" | split -l 25 - testy
ls testy* > running.sh
./parallelize.sh

cat clusters.tsv_* > all_clusters_80_60.tsv
mkdir all_validation_80_60 #See how much of the presence/absence table can be confirmed with a 90_90 search + collapsing paralog hits
cp all_clusters_80_60.tsv validation
cd validation
cut -f1 -d "@" all_clusters_80_60.tsv | awk -F'\t' '{print > $1".tsv"}' #Split into gene-by-gene files
for i in $(ls *tsv | grep -v "cluster"); do sort -nk 2 $i -o $i; done #Sort for comparison

#Now compare with the reported presence/absence pattern:
for i in $(ls *tsv | grep -v "cluster"); do 
    if [ -z "$(comm -3 <(cut -f2 "$i" | sort) <(cut -f2 ../validation_positive/"$i" | sort))" ]; then 
        echo "$i" >> first_pass
    fi
done

cut -f1 -d '.' first_pass | grep -w -F -f - all_clusters_80_60.tsv > validated_clusters_step1.tsv

mkdir paralogs
cp all_clusters_80_60.tsv paralogs
cd paralogs
cut -f1 -d "@" all_clusters_80_60.tsv | sort -u | awk -F'\t' '{print > $1".tsv"}'
for i in $(ls *tsv | grep -v "cluster"); do sort -nk 2 $i -o $i; done
for i in $(ls *tsv | grep -v "cluster"); do 
    if [ -z "$(comm -3 <(cut -f2 "$i" | sort) <(cut -f2 ../../validation_positive/"$i" | sort))" ]; then 
        echo "$i" >> first_pass_paralogs
    fi
done

#Subtract first pass genes from this paralog-containing list

cd ..
cat first_pass paralogs/first_pass_paralogs | sort | uniq -c | grep " 1 " | rev | cut -f1 -d " " | rev | sed "s/\.tsv//g" > paralogs_ids

#Idetify the highest bit-score hit for each paralog gene (only in those target genomes with >1 hit)

for i in $(cat paralogs_ids)
do
echo "grep -w \"$i\" all_clusters_80_60.tsv | cut -f1 -d '@' | sort | uniq -c | awk '(\$1!=1)' | rev | cut -f1 -d ' ' | rev | grep -f - ../500_gffs/500_proteins/ref_vs_500concat.tsv | awk 'BEGIN {max = 0} {if (\$17 > max) {max = \$17; line = \$0}} END {print line}' | cut -f1,2 >> paralog_highest_score.tsv"; done |
awk '{print $0"_"NR}' > running.sh
./parallelize_run.sh

cat paralog_highest_score.tsv_* > paralog_highest_score_hits.tsv && rm paralog_highest_score.tsv_*

#Replace all cases where a paralog has more than one hit against a genome
#With only its best-scoring hit
grep -w -F -f paralogs_ids all_clusters_80_60.tsv > temp
cut -f 1 -d "@" temp | sort | uniq -c | awk '($1!=1)' | rev | cut -f1 -d " " | rev | grep -w -F -vf - temp > validated_clusters_step2.tsv
cat paralog_highest_score_hits.tsv >> validated_clusters_step2.tsv

#Now for the remainder

cat validated_clusters_step1.tsv validated_clusters_step2.tsv | cut -f1 | sort -u | grep -F -w -v -f - all_clusters_80_60.tsv | cut -f1 | sort -u > remaining_genes.txt
cp remaining_genes.txt remainder
cat validated_clusters_step1.tsv validated_clusters_step2.tsv | cut -f1 | sort -u | grep -F -w -v -f - all_clusters_80_60.tsv > remaining_clusters_80_60.tsv
cp remaining_clusters_80_60.txt remainder/
cd remainder/

#paralog-inclusive

cut -f1 -d "@" remaining_clusters_80_60.tsv | sort -u | awk -F'\t' '{print > $1".tsv"}'
#Count number of targets in case of remaining genes
for i in $(ls *tsv | grep -v "remain"); do wc -l $i; done > remainder_target_count.tsv
for i in $(ls *tsv | grep -v "remain"); do wc -l ../../validation_positive/$i; done > remainder_supposed_count.tsv
sort -k2 remainder_target_count.tsv -o remainder_target_count.tsv
sed -i "s/\.\.\/\.\.\/validation_positive\///g" remainder_supposed_count.tsv
sort -k2 remainder_supposed_count.tsv -o remainder_supposed_count.tsv

#Isolate cases where the #strains in which the gene is found is 90% of the #strains in which it's supposed to be found
paste remainder_target_count.tsv remainder_supposed_count.tsv | sed "s/\.tsv//g" | sed "s/ /\t/g" | awk -F '\t' '{OFS=FS}{print $2,$1,$3,$3-$1,$1/$3}' | sort -nrk 5 | awk -F '\t' '($5>0.9)' | cut -f1 > remainder_highquality_ids.txt

#Which one of these has at least one >1 hit against a genome?
grep -w -F -f remainder_highquality_ids.txt remaining_clusters_80_60.tsv | cut -f1 -d "@" | sort | uniq -c | awk '($1!=1)' | rev | cut -f1 -d " " | rev | cut -f1 | sort -u > remainder_paralog_ids.txt
#The rest:
grep -w -F -v -f remainder_paralog_ids.txt remainder_highquality_ids.txt > remainder_nonparalog_ids.txt
grep -w -F -f remainder_nonparalog_ids.txt remaining_clusters_80_60.tsv > validated_clusters_step3.tsv

#Process paralogs:
#Idetify the highest bit-score hit for each paralog gene (only in those target genomes with >1 hit)

for i in $(cat remainder_paralog_ids.txt)
do
echo "grep -w \"$i\" remaining_clusters_80_60.tsv | cut -f1 -d '@' | sort | uniq -c | awk '(\$1!=1)' | rev | cut -f1 -d ' ' | rev | grep -f - ../../500_gffs/500_proteins/ref_vs_500concat.tsv | awk 'BEGIN {max = 0} {if (\$17 > max) {max = \$17; line = \$0}} END {print line}' | cut -f1,2 >> paralog_highest_score.tsv"; done |
awk '{print $0"_"NR}' > running.sh
./parallelize_run.sh

cat paralog_highest_score.tsv_* > paralog_highest_score_hits.tsv && rm paralog_highest_score.tsv_*

#Replace all cases where a paralog has more than one hit against a genome
#With only its best-scoring hit
grep -w -F -f remainder_paralog_ids.txt remaining_clusters_80_60.tsv > temp
cut -f 1 -d "@" temp | sort | uniq -c | awk '($1!=1)' | rev | cut -f1 -d " " | rev | grep -w -F -vf - temp > validated_clusters_step4.tsv
cat paralog_highest_score_hits.tsv >> validated_clusters_step4.tsv

#What the hell let's get the rest

mkdir weird_1648
cat ../validated_clusters_step1.tsv ../validated_clusters_step2.tsv validated_clusters_step3.tsv validated_clusters_step4.tsv | cut -f1 | sort -u | grep -v -w -F -f - ../../30547_total > weird_ids.txt
cp weird_ids.txt weird_1648
#Which one of these has at least one >1 hit against a genome?
cd weird_1648
cp remaining_clusters_80_60.tsv weird_1648
grep -w -F -f weird_ids.txt remaining_clusters_80_60.tsv | cut -f1 -d "@" | sort | uniq -c | awk '($1!=1)' | rev | cut -f1 -d " " | rev | cut -f1 | sort -u > weird_paralog_ids.txt
#The rest:
grep -w -F -v -f weird_paralog_ids.txt weird_ids.txt > weird_nonparalog_ids.txt
grep -w -F -f weird_nonparalog_ids.txt remaining_clusters_80_60.tsv > validated_clusters_step5.tsv

#Process paralogs:
#Idetify the highest bit-score hit for each paralog gene (only in those target genomes with >1 hit)

for i in $(cat weird_paralog_ids.txt)
do
echo "grep -w \"$i\" remaining_clusters_80_60.tsv | cut -f1 -d '@' | sort | uniq -c | awk '(\$1!=1)' | rev | cut -f1 -d ' ' | rev | grep -f - ../../../500_gffs/500_proteins/ref_vs_500concat.tsv | awk 'BEGIN {max = 0} {if (\$17 > max) {max = \$17; line = \$0}} END {print line}' | cut -f1,2 >> paralog_highest_score.tsv"; done |
awk '{print $0"_"NR}' > running.sh
./parallelize_run.sh

cat paralog_highest_score.tsv_* > paralog_highest_score_hits.tsv && rm paralog_highest_score.tsv_*

#Replace all cases where a paralog has more than one hit against a genome
#With only its best-scoring hit
grep -w -F -f weird_paralog_ids.txt remaining_clusters_80_60.tsv > temp
cut -f 1 -d "@" temp | sort | uniq -c | awk '($1!=1)' | rev | cut -f1 -d " " | rev | grep -w -F -vf - temp > validated_clusters_step6.tsv
cat paralog_highest_score_hits.tsv >> validated_clusters_step6.tsv

#Concatenate the four files, cut -f1 -d "@", then reconstruct the p/a table

cd /stor/work/Ochman/hassan/Ecoli_pangenome/all_validation_80_60
cat validated_clusters_step1.tsv | sed "s/$/\tfirstpass/g" > protein_clusters.tsv
cat validated_clusters_step2.tsv | sed "s/$/\tfirstpass_paralogs/g" >> protein_clusters.tsv
cat remainder/validated_clusters_step3.tsv | sed "s/$/\tsecondpass/g" >> protein_clusters.tsv
cat remainder/validated_clusters_step4.tsv | sed "s/$/\tsecondpass_paralogs/g" >> protein_clusters.tsv
cat remainder/weird_1648/validated_clusters_step5.tsv | sed "s/$/\tthirdpass/g" >> protein_clusters.tsv
cat remainder/weird_1648/validated_clusters_step6.tsv | sed "s/$/\tthirdpass_paralogs/g" >> protein_clusters.tsv
sort -u protein_clusters.tsv -o protein_clusters.tsv

#Left: 30547-30216 = 331

#Firstpass: Diamond search 80/60 agrees with presence/absence, no paralog consideration necessary. P/A per original
#Firstpass_paralogs: Diamond search 80/60 agrees with presence/absence, paralog consideration necessary. P/A per original
#Secondpass: Diamond search 80/60 kinda agrees with presence/absence (90% in number of strains), no paralog consideration necessary. P/A per diamond 80/60
#Secondpass_paralogs: Diamond search 80/60 kinda agrees with presence/absence (90% in number of strains), no paralog consideration necessary. P/A per diamond 80/60
#Thirdpass: Everything else that had a 80/60 hit, no paralog consideration necessary. P/A per diamond 80/60
#Fourthpass: Everything else that had a 80/60 hit, paralog consideration necessary. P/A per diamond 80/60
