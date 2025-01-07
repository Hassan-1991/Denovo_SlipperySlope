#Manually assessing alignments:

#First, identify all non-ORFans:

mkdir potential_nonORFans
cd potential_nonORFans
mkdir genus
mkdir species
cd ..
cut -f1,50- ExcelInput.wide.ordered.tsv | grep "candidate" | cut -f1 | sort -u | sed "s/^/cp flanks\//g" | sed "s/$/\*aln potential_nonORFans\/genus\//g" | bash
cd potential_nonORFans/genus

for i in $(ls *aln | cut -f1,2 -d "_")
do
seqkit fx2tab "$i"_mafft.aln | egrep -v "Escherichia|Ecoli@" | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" > "$i"_mafft.reduced.aln
done
rm *mafft.aln

#Detecting species-specific non-ORFans:
cd ../../
cut -f1,45- ExcelInput.wide.ordered.tsv | grep -v "present" | cut -f1 | tail -n+2 | grep -F -f - ExcelInput.wide.ordered.tsv | cut -f1,45- | grep "candidate" | cut -f1 | sort -u | sed "s/^/cp flanks\//g" | sed "s/$/\*aln potential_nonORFans\/species\//g" | bash
cd potential_nonORFans/species/

for i in $(ls *aln | cut -f1,2 -d "_")
do
seqkit fx2tab "$i"_mafft.aln | egrep -v "Ecoli@" | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" > "$i"_mafft.reduced.aln
done
rm *mafft.aln

#After manually detecting non-ORFans, it's time to label the initial list of genus-specific ORFans
#At least those that had at least one pair of flanks in the pangenome

CHJCLABC_02709 - genus_flank_ORFan hobe
PMLIEION_02833 - species_flank_nonORFan hobe
CNLDCNHB_03529 - genus_flank_ORFan hobe
COGJFENE_04297 - genus_flank_ORFan hobe
BALMIOGO_00355 - species_flank_nonORFan hobe
CHJCLABC_00423 - genus_flank_ORFan hobe
OHPFDDCB_04465 - genus_flank_ORFan, species_flank_nonORFan hobe
HJNDJCHN_04802 - species_flank_nonORFan hobe

#CODE
awk -F ',' '{print $2"\t"$1}' interim2 | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/CHJCLABC_02709\tgenus_flank_nonORFan/CHJCLABC_02709\tgenus_flank_ORFan/g" | sed "s/PMLIEION_02833\tgenus_flank_nonORFan, species_flank_ORFan/PMLIEION_02833\tgenus_flank_nonORFan, species_flank_nonORFan/g" | sed "s/CNLDCNHB_03529\tgenus_flank_nonORFan, species_flank_ORFan/CNLDCNHB_03529\tgenus_flank_ORFan, species_flank_ORFan/g" | sed "s/COGJFENE_04297\tgenus_flank_nonORFan, species_flank_ORFan/COGJFENE_04297\tgenus_flank_ORFan, species_flank_ORFan/g" | sed "s/BALMIOGO_00355\tgenus_flank_nonORFan, species_flank_ORFan/BALMIOGO_00355\tgenus_flank_nonORFan, species_flank_nonORFan/g" | sed "s/CHJCLABC_00423\tgenus_flank_nonORFan, species_flank_ORFan/CHJCLABC_00423\tgenus_flank_ORFan, species_flank_ORFan/g" | sed "s/OHPFDDCB_04465\tgenus_flank_nonORFan, species_flank_ORFan/OHPFDDCB_04465\tgenus_flank_ORFan, species_flank_nonORFan/g" | sed "s/HJNDJCHN_04802\tgenus_flank_nonORFan, species_flank_ORFan/HJNDJCHN_04802\tgenus_flank_nonORFan, species_flank_nonORFan/g" | cut -f2 | sort | uniq -c

####ALSO detecting false positives not stuck in conserved syntenic positions####

seqkit fx2tab all_450_CDS.faa | sort -u | sed "s/\t$//g" | sed "s/^/>/g" | grep -F -f genus_specific_ORFans_largeenoughcontigs.txt | sed "s/\t/\n/g" > genus_specific_ORFans_largeenoughcontigs.CDS.faa
cut -f1,45- ExcelInput.wide.ordered.tsv | grep -v "present" | cut -f1 | grep -F -f - ExcelInput.wide.ordered.tsv | tail -n+2 | cut -f1 > species_specific_ORFans_largeenoughcontigs.txt
seqkit fx2tab all_450_CDS.faa | sort -u | sed "s/\t$//g" | sed "s/^/>/g" | grep -F -f species_specific_ORFans_largeenoughcontigs.txt | sed "s/\t/\n/g" > species_specific_ORFans_largeenoughcontigs.CDS.faa

#Genus-specific blasn
blastn -query genus_specific_ORFans_largeenoughcontigs.CDS.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out genus_specific_ORFans_largeenoughcontigs.blastn
#Species-specific blasn
blastn -query species_specific_ORFans_largeenoughcontigs.CDS.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out species_specific_ORFans_largeenoughcontigs.blastn

awk '/^Query=/ {close(file); match($0, /@(.*)\(/, arr); file = arr[1] "_blastn"} {if (file) print > file}' genus_specific_ORFans_largeenoughcontigs.blastn
grep "No hits found" *_blastn | cut -f1,2 -d "_" | sed "s/^/rm /g" | sed "s/$/\*/g" | bash

header=$(head -14 genus_specific_ORFans_largeenoughcontigs.blastn)
footer=$(tail -11 genus_specific_ORFans_largeenoughcontigs.blastn)
for file in *_blastn; do
    # Prepend the header to the file
    { echo "$header"; cat "$file"; } > temp_file && mv temp_file "$file"
    
    # Append the footer to the file
    { cat "$file"; echo "$footer"; } > temp_file && mv temp_file "$file"
done

mkdir genus_specific
mv *_blastn genus_specific
cd genus_specific
export PERL5LIB=/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/lib/
for i in $(ls *_blastn | cut -f1,2 -d "_")
do
echo "/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/bin/mview -in blast ${i}_blastn > ${i}_blastn_mviewed"
done > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh

for i in $(ls *mviewed | cut -f1,2 -d "_")
do
querylength=$(grep -A1 "$i" ../../genus_specific_ORFans_largeenoughcontigs.CDS.faa | grep -v "^>" | awk '{print length($0)}')
ratio=$(tail -n+8 "$i"_blastn_mviewed | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}')
if (( $(echo "$ratio > 0.5" | bc -l) ))
then
tail -n+9 "$i"_blastn_mviewed | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2>50)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear > "$i"_blastn_seq.faa
tail -n+8 "$i"_blastn_mviewed | head -n-3 | sed "1s/^/g /g" | awk '{print $2,$NF}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | head -2 | seqkit fx2tab | cut -f2- -d "@" | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" >> "$i"_blastn_seq.faa
fi
done

#Reduce redundancy:

egrep -v "Escherichia|@" ../../all_contig_protein_taxonomy* > genus_contig_protein_taxonomy.tsv
for i in $(ls *_blastn_seq.faa | cut -f1,2 -d "_")
do
grep "^>" "$i"_blastn_seq.faa | tr -d ">" | head -n-1 | sort -u | grep -w -F -f - genus_contig_protein_taxonomy.tsv | sort -k1 > interim
seqkit fx2tab "$i"_blastn_seq.faa | sort -k1 | join -1 1 -2 1 - interim | awk '{print $1,$2":"$3}' | awk '!seen[$2]++' | sed "s/:/ /g" | awk '{print $3":"$1"\t"$2}' | sed "s/^/>/g" | sed "s/\t/\n/g" > "$i"_mafft_input.faa
tail -2 "$i"_blastn_seq.faa >> "$i"_mafft_input.faa
grep -A1 "$i" ../../genus_specific_ORFans_largeenoughcontigs.CDS.faa | seqkit fx2tab | sed "s/\t$//g" | cut -f2- -d "@" | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/>/>FULL_/g" | sed "s/\t/\n/g" >> "$i"_mafft_input.faa
done

for i in $(ls *_mafft_input.faa | cut -f1,2 -d "_")
do
mafft --auto "$i"_mafft_input.faa > "$i"_mafft.aln
done

#SPECIES#
mkdir species_specific
mv species_specific_ORFans_largeenoughcontigs.blastn species_specific
awk '/^Query=/ {close(file); match($0, /@(.*)\(/, arr); file = arr[1] "_blastn"} {if (file) print > file}' species_specific_ORFans_largeenoughcontigs.blastn
grep "No hits found" *_blastn | cut -f1,2 -d "_" | sed "s/^/rm /g" | sed "s/$/\*/g" | bash

header=$(head -14 species_specific_ORFans_largeenoughcontigs.blastn)
footer=$(tail -11 species_specific_ORFans_largeenoughcontigs.blastn)
for file in *_blastn; do
    # Prepend the header to the file
    { echo "$header"; cat "$file"; } > temp_file && mv temp_file "$file"
    
    # Append the footer to the file
    { cat "$file"; echo "$footer"; } > temp_file && mv temp_file "$file"
done

for i in $(ls *_blastn | cut -f1,2 -d "_")
do
echo "/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/bin/mview -in blast ${i}_blastn > ${i}_blastn_mviewed"
done > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh

for i in $(ls *mviewed | cut -f1,2 -d "_")
do
querylength=$(grep -A1 "$i" ../../genus_specific_ORFans_largeenoughcontigs.CDS.faa | grep -v "^>" | awk '{print length($0)}')
ratio=$(tail -n+8 "$i"_blastn_mviewed | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}')
if (( $(echo "$ratio > 0.5" | bc -l) ))
then
tail -n+9 "$i"_blastn_mviewed | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2>50)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear > "$i"_blastn_seq.faa
tail -n+8 "$i"_blastn_mviewed | head -n-3 | sed "1s/^/g /g" | awk '{print $2,$NF}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | head -2 | seqkit fx2tab | cut -f2- -d "@" | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" >> "$i"_blastn_seq.faa
fi
done

#Reduce redundancy:

egrep "Escherichia" ../../all_contig_protein_taxonomy* > species_contig_protein_taxonomy.tsv

for i in $(ls *_blastn_seq.faa | cut -f1,2 -d "_")
do
grep "^>" "$i"_blastn_seq.faa | tr -d ">" | head -n-1 | sort -u | grep -w -F -f - species_contig_protein_taxonomy.tsv | sort -k1 > interim
seqkit fx2tab "$i"_blastn_seq.faa | sort -k1 | join -1 1 -2 1 - interim | awk '{print $1,$2":"$3}' | awk '!seen[$2]++' | sed "s/:/ /g" | awk '{print $3":"$1"\t"$2}' | sed "s/^/>/g" | sed "s/\t/\n/g" > "$i"_mafft_input.faa
tail -2 "$i"_blastn_seq.faa >> "$i"_mafft_input.faa
grep -A1 "$i" ../../genus_specific_ORFans_largeenoughcontigs.CDS.faa | seqkit fx2tab | sed "s/\t$//g" | cut -f2- -d "@" | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/>/>FULL_/g" | sed "s/\t/\n/g" >> "$i"_mafft_input.faa
done

for i in $(ls *_mafft_input.faa | cut -f1,2 -d "_")
do
mafft --auto "$i"_mafft_input.faa > "$i"_mafft.aln
done

######IN PARALLEL######

#Make a better getorf database:
getorf -sequence Escherichia_excluded.genomes.faa -outseq Escherichia_excluded.genomes.getorf.bact -table 11 -minsize 30 -find 3


#How to settle cases of missing front and back halves?

#Extract 300bp from both edges from target sequences:

for i in $(ls *mafft.aln | cut -f1,2 -d "_" | grep "BECHDPHP_04018")
do
#no /rc, +
seqkit fx2tab "$i"_mafft.aln | head -n-2 | cut -f2 -d ":" | sort -u | sed "s/$/:/g" |
grep -F -f - /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/backup_2/flanks/"$i"_blastn_mviewed |
awk '($2!~"/rc"&&$7=="+")' | awk '{print $2,$11}' |
sed "s/:/\t/g" | rev | sed "s/-/\t/" | rev |
awk '{print "samtools faidx /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/" $1,$1":"$2+$4-1-100"-"$2+$5-1+100}' | bash > "$i".interim
#no /rc, -
seqkit fx2tab "$i"_mafft.aln | head -n-2 | cut -f2 -d ":" | sort -u | sed "s/$/:/g" |
grep -F -f - /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/backup_2/flanks/"$i"_blastn_mviewed |
awk '($2!~"/rc"&&$7=="-")' | awk '{print $2,$11}' |
sed "s/:/\t/g" | rev | sed "s/-/\t/" | rev |
awk '{print "samtools faidx -i /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/" $1,$1":"$2+$5-1-100"-"$2+$4-1+100}' | bash >> "$i".interim
#/rc, +
seqkit fx2tab "$i"_mafft.aln | head -n-2 | cut -f2 -d ":" | sort -u | sed "s/$/:/g" |
grep -F -f - /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/backup_2/flanks/"$i"_blastn_mviewed |
awk '($2~"/rc"&&$7=="+")' | awk '{print $2,$11}' |
sed "s/:/\t/g" | rev | sed "s/-/\t/" | rev | sed "s/\/rc//g" |
awk '{print "samtools faidx -i /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/" $1,$1":"$3-$5+1-100"-"$3-$4+1+100}' | bash >> "$i".interim
#/rc, -
seqkit fx2tab "$i"_mafft.aln | head -n-2 | cut -f2 -d ":" | sort -u | sed "s/$/:/g" |
grep -F -f - /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/backup_2/flanks/"$i"_blastn_mviewed |
awk '($2~"/rc"&&$7=="-")' | awk '{print $2,$11}' |
sed "s/:/\t/g" | rev | sed "s/-/\t/" | rev | sed "s/\/rc//g" |
awk '{print "samtools faidx /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/" $1,$1":"$3-$4+1-100"-"$3-$5+1+100}' | bash >> "$i".interim
#extract contig names and taxonomy
seqkit fx2tab "$i"_mafft.aln | head -n-2 | cut -f1,2 -d ":" | sed "s/:/\t/g" | sort -u | sort -k2 > "$i".interim.2
seqkit fx2tab "$i".interim | sed "s/:/\t/" | sort -k1 | join -1 1 -2 2 - "$i".interim.2 | awk '{print ">"$NF":"$1":"$2"\t"$3}' | sed "s/\t$//g" | sed "s/\t/\n/g" > "$i"_secondstep.faa
seqkit fx2tab "$i"_secondstep.faa | sed "s/:/\t/" | awk -F '\t' '{print $1"%"$3,$2}' | awk '!seen[$1]++' | sed "s/%/\t/g" | awk '{print ">"$1":"$3,$2}' | sed "s/ /\n/g" > "$i"_secondstep.nr.faa
seqkit fx2tab "$i"_mafft.aln | grep "FULL" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" >> "$i"_secondstep.nr.faa
mafft --auto "$i"_secondstep.nr.faa > "$i"_secondstep.aln
done

#no /rc, +

seqkit fx2tab ACOBLJOG_01406_mafft.aln | head -n-2 | cut -f2 -d ":" | sort -u | sed "s/$/:/g" |
grep -F -f - /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/backup_2/flanks/ACOBLJOG_01406_blastn_mviewed |
awk '($2!~"/rc"&&$7=="+")' | awk '{print $2,$11}' |
sed "s/:/\t/g" | rev | sed "s/-/\t/" | rev |
awk '{print "samtools faidx /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/" $1,$1":"$2+$4-1-300"-"$2+$5-1+300}' | bash > ACOBLJOG_01406.interim

#no /rc, -

seqkit fx2tab ACOBLJOG_01406_mafft.aln | head -n-2 | cut -f2 -d ":" | sort -u | sed "s/$/:/g" |
grep -F -f - /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/backup_2/flanks/ACOBLJOG_01406_blastn_mviewed |
awk '($2!~"/rc"&&$7=="-")' | awk '{print $2,$11}' |
sed "s/:/\t/g" | rev | sed "s/-/\t/" | rev |
awk '{print "samtools faidx -i /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/" $1,$1":"$2+$5-1-300"-"$2+$4-1+300}' | bash >> ACOBLJOG_01406.interim

#/rc, +

seqkit fx2tab ACOBLJOG_01406_mafft.aln | head -n-2 | cut -f2 -d ":" | sort -u | sed "s/$/:/g" |
grep -F -f - /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/backup_2/flanks/ACOBLJOG_01406_blastn_mviewed |
awk '($2~"/rc"&&$7=="+")' | awk '{print $2,$11}' |
sed "s/:/\t/g" | rev | sed "s/-/\t/" | rev | sed "s/\/rc//g" |
awk '{print "samtools faidx -i /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/" $1,$1":"$3-$5+1-300"-"$3-$4+1+300}' | bash >> ACOBLJOG_01406.interim

#/rc, -

seqkit fx2tab ACOBLJOG_01406_mafft.aln | head -n-2 | cut -f2 -d ":" | sort -u | sed "s/$/:/g" |
grep -F -f - /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/backup_2/flanks/ACOBLJOG_01406_blastn_mviewed |
awk '($2~"/rc"&&$7=="-")' | awk '{print $2,$11}' |
sed "s/:/\t/g" | rev | sed "s/-/\t/" | rev | sed "s/\/rc//g" |
awk '{print "samtools faidx /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/" $1,$1":"$3-$4+1-300"-"$3-$5+1+300}' | bash >> ACOBLJOG_01406.interim

#extract contig names and taxonomy
seqkit fx2tab ACOBLJOG_01406_mafft.aln | head -n-2 | cut -f1,2 -d ":" | sed "s/:/\t/g" | sort -u | sort -k2 > ACOBLJOG_01406.interim.2

#
seqkit fx2tab ACOBLJOG_01406.interim | sed "s/:/\t/" | sort -k1 | join -1 1 -2 2 - ACOBLJOG_01406.interim.2 | awk '{print ">"$NF":"$1":"$2"\t"$3}' | sed "s/\t$//g" | sed "s/\t/\n/g" > ACOBLJOG_01406_secondstep.faa
seqkit fx2tab ACOBLJOG_01406_secondstep.faa | sed "s/:/\t/" | awk -F '\t' '{print $1"%"$3,$2}' | awk '!seen[$1]++' | sed "s/%/\t/g" | awk '{print ">"$1":"$3,$2}' | sed "s/ /\n/g" > ACOBLJOG_01406_secondstep.nr.faa
seqkit fx2tab ACOBLJOG_01406_mafft.aln | grep "FULL" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" >> ACOBLJOG_01406_secondstep.nr.faa

mafft --auto ACOBLJOG_01406_secondstep.nr.faa > ACOBLJOG_01406_secondstep.aln

#####

#Do a blastn search with all ORFans

seqkit fx2tab all_450_CDS.faa | sort -u | sed "s/\t$//g" | sed "s/^/>/g" | grep -F -f genus_specific_ORFans_largeenoughcontigs.txt | sed "s/\t/\n/g" > genus_specific_ORFans_largeenoughcontigs.CDS.faa

blastn -query genus_specific_ORFans_largeenoughcontigs.CDS.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out geneflanks_extragenus_blastn

blastn -query genus_specific_ORFans_largeenoughcontigs.CDS.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out genus_specific_ORFans_largeenoughcontigs.smallwordblastn -word_size 7

