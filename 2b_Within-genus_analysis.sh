mkdir noncoli_Escherichia_flankanalysis_updated
cp Ecoli_step4_ORFans.final.proxflanks.faa noncoli_Escherichia_flankanalysis_updated
cp Ecoli_step4_ORFans.final.geneflanks.faa noncoli_Escherichia_flankanalysis_updated

#1545 genes have proximal flanks (500bp)
#1469 has gene flanks
#Altogether, 1545 genes have at least some kinda flank

blastn -query Ecoli_step4_ORFans.final.proxflanks.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Ecoli_step4_ORFans_proxflanks_blastn
cat Ecoli_step4_ORFans_proxflanks_blastn | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | cut -f 1 -d "@" | awk -F '_' '{OFS="\t"}{print $1,$0}' | sed "s/left500_//g" | sed "s/right500_//g" | awk -F '\t' '{OFS=""}{print $1,"\t",$2,"%",$3}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep "left.*right" | cut -f 1 | awk -F'%' '{ values[$1] = (values[$1] == "" ? $2 : values[$1] ", " $2) } END { for (value in values) { print value "\t" values[value] } }' > Ecoli_step4_ORFans_proxflanks_targets.txt

blastn -query Ecoli_step4_ORFans.final.geneflanks.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Ecoli_step4_ORFans_geneflanks_blastn
cat Ecoli_step4_ORFans_geneflanks_blastn | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | cut -f 1 -d "@" | sed "s/_up_/%up%/g" | sed "s/_down_/%down%/g" | awk -F "\t" '{OFS=""}{print $1,"%",$2}' | awk -F '%' '{OFS=""}{print $2,"_",$3,"\t",$1,"%",$4}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | egrep "up.*down|down.*up" | cut -f 1 | awk -F'%' '{ values[$1] = (values[$1] == "" ? $2 : values[$1] ", " $2) } END { for (value in values) { print value "\t" values[value] } }' > Ecoli_step4_ORFans_geneflanks_targets.txt

for i in $(cut -f 1 -d "(" Ecoli_step4_ORFans_proxflanks_targets.txt)
do
grep "$i(" Ecoli_step4_ORFans_proxflanks_targets.txt | sed "s/,/\n/g" | sed "s/)\t/)\n/g" | sed "s/^ *//g" | tail -n+2 > "$i"_Ecoli_step4_ORFans_proxflanks_targetlist.txt
grep -f "$i"_Ecoli_step4_ORFans_proxflanks_targetlist.txt Ecoli_step4_ORFans_proxflanks_blastn | grep "$i(" | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | awk -F'\t' 'BEGIN { OFS = "\t" } { sub(/@.*$/, "", $2); print }' | cut -f 2- -d "_" | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > "$i"_Ecoli_step4_ORFans_proxflanks_intervalinfo
done

for i in $(cut -f 1 Ecoli_step4_ORFans_geneflanks_targets.txt)
do
grep -P "$i\t" Ecoli_step4_ORFans_geneflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 > "$i"_Ecoli_step4_ORFans_geneflanks_targetlist.txt
grep -f "$i"_Ecoli_step4_ORFans_geneflanks_targetlist.txt Ecoli_step4_ORFans_geneflanks_blastn | grep "$i"_ | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | awk -F'\t' 'BEGIN { OFS = "\t" } { sub(/@.*$/, "", $2); print }' | sed "s/_up/%up/g" | sed "s/_down/%down/g" | cut -f 2- -d "%" | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}'  | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > "$i"_Ecoli_step4_ORFans_geneflanks_intervalinfo
done

for i in $(ls *prox*info | cut -f1,2 -d "_")
do
    value2=$(grep "$i" ../Ecoli_step4_ORFans.final.besthits.gtf | awk -F '\t' '{print $5-$4+1}')
    sed -i "s/$/ $value2/" "$i"_Ecoli_species_specific_ORFans_step3_proxflanks_intervalinfo
done

#sidebar: meaningful flanks for comparative analysis between ORFans and non-ORFans
for i in $(ls *prox*info | cut -f1,2 -d "_"); do awk '(($4<($6+600))&&($4>($6*0.5)))' "$i"_Ecoli_species_specific_ORFans_step3_proxflanks_intervalinfo > "$i"_Ecoli_species_specific_ORFans_step3_proxflanks_intervalinfo_2; done
for i in $(ls *gene*info | cut -f1,2 -d "_"); do awk '(($4<($6+2000))&&($4>($6*0.5)))' "$i"_Ecoli_species_specific_ORFans_step3_geneflanks_intervalinfo > "$i"_Ecoli_species_specific_ORFans_step3_geneflanks_intervalinfo_2; done
find . -type f -empty -delete

#At this point, the two streams can hopefully join

for i in $(ls *intervalinfo_2 | cut -f1,2 -d "_" | sort -u); do ls *"$i"_*intervalinfo_2 | sed "s/^/cat /g" | bash >> "$i"_joint_intervalcoords; done
