#For any de novo gene candidate, concatenate all of its blastn hit sequences into one file, then run macse on them.
#Do this separately for tblastn, fwiw

for i in $(cat 21_highest_prob_candidates.txt)
do
grep -A1 "$i" ../step2_ORFans_wflanks_besthits.faa > "$i"_all_blastn_seq.faa
head -n-4 "$i"_blastn_seq.faa >> "$i"_all_blastn_seq.faa
sed "s/>/>fergusonii_/g" ../flanks_species/"$i"*fergusonii_blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
sed "s/>/>ruysiae_/g" ../flanks_species/"$i"*ruysiae_blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
sed "s/>/>marmotae_/g" ../flanks_species/"$i"*marmotae_blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
sed "s/>/>whittamii_/g" ../flanks_species/"$i"*whittamii_blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
sed "s/>/>albertii_/g" ../flanks_species/"$i"*albertii_blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
sed "s/>/>outside_/g" ../flanks_genus/"$i"*blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
mafft --auto "$i"_all_blastn_seq.faa > "$i"_all_blastn_seq.aln
done

for i in $(cut -f1,2,4,6,8,10,12 interim_tallying | sed "s/$/\t/g" | sed "s/no flanks/no flank/g" | sed "s/\tn\t/\tpresent\t/g" | awk -F '\t' '($2!="present")' | awk -F '\t' '{OFS="\t";for (i=2; i<=7; i++) {if ($i != "untraceable" && $i != "no flank" && $i != "present") {$i = "non-coding";}} print $0}' | cut -f2- | sort | uniq -c | sort -nrk1 | sed "s/^ *//g" | sed "s/ /\t/" | cat -n | sed "s/^ *//g" | cut -f3- | egrep "trace|coding" | grep -f - interim_tallying_3 | cut -f1)
do
grep -A1 "$i" ../step2_ORFans_wflanks_besthits.faa > "$i"_all_blastn_seq.faa
head -n-4 "$i"_blastn_seq.faa >> "$i"_all_blastn_seq.faa
sed "s/>/>fergusonii_/g" ../flanks_species/"$i"*fergusonii_blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
sed "s/>/>ruysiae_/g" ../flanks_species/"$i"*ruysiae_blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
sed "s/>/>marmotae_/g" ../flanks_species/"$i"*marmotae_blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
sed "s/>/>whittamii_/g" ../flanks_species/"$i"*whittamii_blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
sed "s/>/>albertii_/g" ../flanks_species/"$i"*albertii_blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
sed "s/>/>outside_/g" ../flanks_genus/"$i"*blastn_seq.faa | head -n-4 >> "$i"_all_blastn_seq.faa
mafft --auto "$i"_all_blastn_seq.faa > "$i"_all_blastn_seq.aln
done




#ls AELFKIGO_03201_all_blastn_seq.faa | awk '{OFS=""}{print "java -jar \/stor\/work\/Ochman\/hassan\/tools\/macse_v2.06.jar -prog alignSequences -fs 1000 -fs_term 1000 -stop 100 -seq ",$1," -out_NT ",$1,"_NT.aln -out_AA ",$1,"_AA.aln"}' | bash
