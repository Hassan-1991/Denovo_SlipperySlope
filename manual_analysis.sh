#For any de novo gene candidate, concatenate all of its blastn hit sequences into one file, then run mafft on them.

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

#Observe these visually with SeaView to identify potential de novo genes
#42 identified from the first 126 surveyed
#Identify their point of origin in the pangenome
#For this, make a presence/absence table for each 419 strains
#The presence/absence for blastp is already there (merged_output_3.tsv)
#Make presence absence for non-coding

#Since de novo calls are made based on visual observation of alignments, mark their status in the NR set of identifiers

8939_4#2	non-coding
SRR7265856_trimmed	non-coding
SRR7187794_trimmed	non-coding
JSHJ00000000	non-coding
JSMU00000000	non-coding
11658_3#10	non-coding
11791_7#51	present
esc_ga8815aa_as	non-coding
esc_pa9513aa_as	non-coding

sort -k1 AABIHECD_00503_nr_identifiers.tsv -o AABIHECD_00503_nr_identifiers.tsv

#Now generate a cluster file of collapsed identifiers
head -n-4 AABIHECD_00503_blastn_seq.faa_2 | seqkit rmdup -s - > AABIHECD_00503_blastn_seq.faa_dedup.faa
python3 seqkit_duplicate_parser.py AABIHECD_00503_blastn_seq.faa_2 AABIHECD_00503_blastn_seq.faa_dedup.faa | sed "s/\t/@/g" | cut -f1,3 -d "@" | sed "s/@/\t/g" | sort -u | sort -k1 | join -1 1 -2 1 - AABIHECD_00503_nr_identifiers.tsv | cut -f2- -d " " | sed "s/ /\t/g" > AABIHECD_00503_419_presence_absence.tsv
cut -f1 AABIHECD_00503_419_presence_absence.tsv | grep -w -vf - ../../../../419_isolate_names.txt | grep -v "Strain" | sed "s/$/\tabsent/g" >> AABIHECD_00503_419_presence_absence.tsv
sort -k1 AABIHECD_00503_419_presence_absence.tsv -o AABIHECD_00503_419_presence_absence.tsv
#Lineage assignment
sort -k1 ../../../../419_ipp_lineagedesignations.tsv | cut -f1,2,5 | join -1 1 -2 1 - AABIHECD_00503_419_presence_absence.tsv | cut -f3,4 -d " " | sort | uniq -c | sed "s/^ *//g" | awk '{print $2,$1,$3}' | sort -nk 1 | sed "s/ /,/g" > AABIHECD_00503_42_lineage_presence_absence.csv
python3 lineage_presenceabsence_parser.py AABIHECD_00503_42_lineage_presence_absence.csv | awk '{print $2,$3,$4,$5}' | sed "s/ /,/g" > AABIHECD_00503_42_lineage_presence_absence.fraction.csv

#Newick file of all lineages:

A,B1,E:

(((((((((((((((((ESC GA9113AA AS,ESC JA4591AA AS),ESC PA2004AA AS),SRR7283764 TRIMMED),(ESC PA2318AA AS,SRR7367225 TRIMMED)),SRR5031203 TRIMMED),((SRR7186949 TRIMMED,ESC RA4296AA AS),7748 7#47)),((7790 1#64,ESC HA6177AA AS),ESC QA2206AA AS)),SRR7186910 TRIMMED),(ESC HA3950AA AS,SRR7215958 TRIMMED)),11791 3#6),SRR7290889 TRIMMED),ESC CA5569AA AS),SRR5006334 TRIMMED),ESC GA9425AA AS),((SRR7236613 TRIMMED,ESC HA3288AA AS),ESC JA4735AA AS)),((ESC QA2893AA AS,SRR7367461 TRIMMED),ESC JA5404AA AS)),((((((11657 6#33,ESC HA7438AA AS),11657 5#60),(((((((((11679 5#79,11679 7#17),ESC JA4720AA AS),11679 7#36),JSNL00000000),12045 3#57),((11679 6#14,ESC HA8660AA AS),11679 7#78)),ESC PA9513AA AS),((11658 3#8,ESC BA1344AA AS),11657 5#57))))),(((SRR7251057 TRIMMED,ESC RA2561AA AS),ESC BA5409AA AS),(59235 G01 CONTIGS HGAP 4 0,ESC HA8654AA AS))))

B2,D,F:
((((((11657 6#33,ESC HA7438AA AS),11657 5#60),(((((((((11679 5#79,11679 7#17),ESC JA4720AA AS),11679 7#36),JSNL00000000),12045 3#57),((11679 6#14,ESC HA8660AA AS),11679 7#78)),ESC PA9513AA AS),((11658 3#8,ESC BA1344AA AS),11657 5#57)))),(((SRR7251057 TRIMMED,ESC RA2561AA AS),ESC BA5409AA AS),(59235 G01 CONTIGS HGAP 4 0,ESC HA8654AA AS))))
D:

B2:



