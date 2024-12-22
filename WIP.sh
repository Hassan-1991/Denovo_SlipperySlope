#Manually assessing alignments:

#First, get rid of all non-ORFans:

mkdir potential_nonORFans
cut -f1,50- ExcelInput.wide.ordered.tsv | grep "candidate" | cut -f1 | sort -u | sed "s/^/cp flanks\//g" | sed "s/$/\*aln potential_nonORFans/g" | bash
cd potential_nonORFans

for i in $(ls *aln | cut -f1,2 -d "_")
do
seqkit fx2tab "$i"_mafft.aln | egrep -v "Escherichia|Ecoli@" | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" > "$i"_mafft.reduced.aln
done


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

