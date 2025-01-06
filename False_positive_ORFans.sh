#Paste genus, without synteny analysis results in: fp_withoutsynteny_genus_raw

grep "non-ORFan" fp_withoutsynteny_genus_raw | cut -f1,2 -d "_" | sort -u | sed "s/$/,genus_flank_nonORFan_wo/g" > false_positive.interim
grep -v "non-ORFan" fp_withoutsynteny_genus_raw | cut -f1,2 -d "_" | sort -u | sed "s/$/,genus_flank_ORFan_wo/g" >> false_positive.interim

#Paste species, without synteny analysis results in: fp_withoutsynteny_species_raw

grep "non-ORFan" fp_withoutsynteny_species_raw | cut -f1,2 -d "_" | sort -u | sed "s/$/,species_flank_nonORFan_wo/g" >> false_positive.interim
grep -v "non-ORFan" fp_withoutsynteny_species_raw | cut -f1,2 -d "_" | sort -u | sed "s/$/,species_flank_ORFan_wo/g" >> false_positive.interim

#Paste genus, with synteny analysis results in: fp_withsynteny_genus_raw

grep "non-ORFan" fp_withsynteny_genus_raw | cut -f1,2 -d "_" | sort -u | sed "s/$/,genus_flank_nonORFan_with/g" >> false_positive.interim
grep -v "non-ORFan" fp_withsynteny_genus_raw | cut -f1,2 -d "_" | sort -u | sed "s/$/,genus_flank_ORFan_with/g" >> false_positive.interim

#Paste genus, with synteny analysis results in: fp_withsynteny_species_raw

grep "non-ORFan" fp_withsynteny_species_raw | cut -f1,2 -d "_" | sort -u | sed "s/$/,species_flank_nonORFan_with/g" >> false_positive.interim
grep -v "non-ORFan" fp_withsynteny_species_raw | cut -f1,2 -d "_" | sort -u | sed "s/$/,species_flank_ORFan_with/g" >> false_positive.interim

#Makeshift R plot fodder:
grep "nonORFan" false_positive.interim | grep "genus" | awk -F ',' '{print $2"\t"$1}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | cut -f2 | sort | uniq -c | sort -nrk1
grep "nonORFan" false_positive.interim | grep "species" | awk -F ',' '{print $2"\t"$1}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | cut -f2 | sort | uniq -c | sort -nrk1
