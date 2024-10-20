split -l2 step2_ORFans_wflanks_besthits.faa
for i in x*; do rename=$(grep "^>" $i | cut -f1 -d "(" | tr -d ">"); mv $i $rename; done

for i in $(grep "^>" step2_ORFans_wflanks_besthits.faa | cut -f1 -d "(" | tr -d ">")
do
echo "blastn -query ${i} -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out ${i}_blastn -word_size 7"
done #Run in parallel

#Convert blast results to mviewed format

export PERL5LIB=/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/lib/
for i in $(grep "^>" step2_ORFans_wflanks_besthits.faa | cut -f1 -d "(" | tr -d ">")
do
echo "/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/bin/mview -in blast ${i}_blastn > ${i}_mview"
done #Run in parallel

#First, remove all cases where there's at least one identical protein match
#The for loop identifies cases with 100% qcov

for i in $(wc -l *mview | grep -v "total" | grep -v " 1 " | rev | cut -f1 -d " " | rev | sed "s/_mview//g")
do
querylength=$(cat $i | grep -v "^>" | awk '{print length($0)}')
ratio=$(tail -n+8 "$i"_mview | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}')
if (( $(echo "$ratio == 1" | bc -l) ))
then
tail -n+9 "$i"_mview | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2==100)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | seqkit rmdup -s | linear > "$i"_exclude_interim
fi
done

#Remove all files with no perfect qcov hits
wc -l *_exclude_interim | grep -v "total" | grep " 0 " | rev | cut -f1 -d " " | rev | sed "s/^/rm /g" | bash

#grep the protein sequence against its targets, in case it matches - exclude
for i in $(ls *_exclude_interim | rev | cut -f 3- -d "_" | rev)
do
/stor/work/Ochman/hassan/tools/faTrans -stop "$i"_exclude_interim "$i"_exclude_interim.prot
/stor/work/Ochman/hassan/tools/faTrans -stop "$i" "$i".prot
linear *_exclude_interim.prot
exactmatch_check=$(grep -v "^>" "$i".prot | tr -d "\n" | sed 's/$/\n/' | grep -f - "$i"_exclude_interim.prot | sort -u | wc -l)
echo "$i,$exactmatch_check" >> exactmatch_candidates.csv
done

awk -F ',' '($2>0)' exactmatch_candidates.csv | cut -f1 -d ',' | sort -u | sed "s/$/_/g" > blastn_tobeexcluded

#For the remainder with decent 50/70 alignment, make macse

for i in $(wc -l *mview | grep -v "total" | grep -v " 1 " | rev | cut -f1 -d " " | rev | sed "s/_mview//g")
do
querylength=$(cat $i | grep -v "^>" | awk '{print length($0)}')
ratio=$(tail -n+8 "$i"_mview | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}')
if (( $(echo "$ratio > 0.5" | bc -l) ))
then
tail -n+9 "$i"_mview | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2>70)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | seqkit rmdup -s | linear > "$i"_blastn_seq.faa
tail -n+8 "$i"_mview | head -n-3 | sed "1s/^/g /g" | awk '{print $2,$NF}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | head -2 >> "$i"_blastn_seq.faa
cat $i | sed "s/>/>FULL_/g" >> "$i"_blastn_seq.faa
fi
done

ls *blastn_seq.faa | grep -vf blastn_tobeexcluded | sed "s/^/wc -l /g" | bash | awk '($1!=4)' | rev | cut -f1 -d " " | rev | awk '{OFS=""}{print "java -jar \/stor\/work\/Ochman\/hassan\/tools\/macse_v2.06.jar -prog alignSequences -fs 1000 -fs_term 1000 -stop 100 -seq ",$1," -out_NT ",$1,"_NT.aln -out_AA ",$1,"_AA.aln"}' > running.sh

##########
##########
##########
##########

#tfasty:

