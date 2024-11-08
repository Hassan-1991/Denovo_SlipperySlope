#Starting raw materials:
#1. gff files
#2. genomes extracted from them
#3. gtf files extracted from them
#4. Lineage designations from paper supp - 500_ipp_lineagedesignations.tsv

#All operations in /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline

#All 450 this time

mkdir /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/450_cds
mkdir /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/450_proteins
for i in $(awk -F '\t' '($5<45)' /stor/work/Ochman/hassan/Ecoli_pangenome/500_ipp_lineagedesignations.tsv | cut -f1 | sort -u)
do
awk -F '\t' '($3=="CDS")' /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/"$i".gtf | gtf2bed | bedtools getfasta -s -name -fi /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_genomes/"$i".fasta -bed - | sed "s/>/>"$i"@/g" > /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/450_cds/"$i"_cds.faa
/stor/work/Ochman/hassan/tools/faTrans -stop /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/450_cds/"$i"_cds.faa /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/450_proteins/"$i".prot.faa
done

#Do clustering solely based on mmseqs2 criteria

cat /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/450_proteins/*prot.faa > all_450_proteins.faa
cat /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/450_cds/*.faa > all_450_CDS.faa
linear all_450_CDS.faa
mmseqs createdb all_450_proteins.faa all_450_proteins
mmseqs search all_450_proteins all_450_proteins resultDB tmp --min-seq-id 0.8 -c 0.7 --cov-mode 2
mmseqs convertalis all_450_proteins all_450_proteins resultDB resultDB.m8
mmseqs linclust all_450_proteins clusterDB tmp --min-seq-id 0.8 -c 0.7 --cov-mode 2
mmseqs createtsv all_450_proteins all_450_proteins clusterDB all_450_proteins.clusters.tsv

#Get representative protein sequences:
cut -f1 all_450_proteins.clusters.tsv | sort -u > clustered_ids.txt
faSomeRecords all_450_proteins.faa clustered_ids.txt clustered_representative_proteins.faa

#I don't have ORF databases for all 500 genomes
awk -F '\t' '($5<45)' /stor/work/Ochman/hassan/Ecoli_pangenome/500_ipp_lineagedesignations.tsv | cut -f1 | sort -u | sed "s/^/cat \/stor\/work\/Ochman\/hassan\/Ecoli_pangenome\/500_gffs\/500_genomes\//g"  | sed "s/$/.fasta >> all_450_genomes.faa/g" | bash
getorf -sequence all_450_genomes.faa -outseq all_450_genomes.getorf.all -table 1 -minsize 30 -find 3
seqkit fx2tab all_450_genomes.getorf.all | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/" > all_450_genomes.getorf.ATG_TTG_GTG
/stor/work/Ochman/hassan/tools/faTrans -stop all_450_genomes.getorf.ATG_TTG_GTG all_450_genomes.getorf.ATG_TTG_GTG.prot.faa
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in all_450_genomes.getorf.ATG_TTG_GTG.prot.faa --db all_450_genomes.getorf.ATG_TTG_GTG.prot

#Likewise for annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in all_450_proteins.faa --db all_450_proteins

#Automatic searches:

#Outside genus, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q clustered_representative_proteins.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_proteins/04072024_nonEscherichia_protein --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_GBRS_annotated.tsv -k 0 -b8 -c1
#Outside genus, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q clustered_representative_proteins.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes.getorf.noCTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_GBRS_ORFs.tsv -k 0 -b8 -c1
#Outside species, annotated 
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q clustered_representative_proteins.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_annotated.tsv -k 0 -b8 -c1
#Outside species, ORFs_ATG
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q clustered_representative_proteins.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes.getorf.ATG --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_ORFs_ATG.tsv -k 0 -b8 -c1
#Outside species, ORFs_TTG_GTG (running)
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q clustered_representative_proteins.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes.getorf.TTG_GTG --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_ORFs_TTG_GTG.tsv -k 0 -b8 -c1
#Across pangenome, annotated (running)
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q clustered_representative_proteins.faa -d /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_annotated.tsv -k 0 -b8 -c1
#Across pangenome, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q clustered_representative_proteins.faa -d /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_ORFs.tsv -k 0 -b8 -c1

#Restrict analysis to genus-specific ORFans:
cat all_proteins_vs_GBRS_annotated.tsv all_proteins_vs_GBRS_ORFs.tsv | awk -F '\t' '($5>60&&$16<0.001)' | cut -f1 | sort -u > step1_genusspecific_nonORFan.txt
grep -vf step1_genusspecific_nonORFan.txt clustered_representative_proteins.faa | grep "^>" | tr -d ">" > step1_genusspecific_ORFan.txt
faSomeRecords clustered_representative_proteins.faa step1_genusspecific_ORFan.txt step1_genusspecific_ORFans.faa

#Get flanks for these proteins

#I need to know whether the upstream and downstream flanks of these proteins are consistent, i.e. belong to the same cluster:

mkdir flank_tangent
cp all_450_proteins.clusters.tsv flank_tangent
cp step1_genusspecific_ORFan.txt flank_tangent
cp all_450_CDS.faa flank_tangent
#Make a list of each cluster to which upstream and downstream flanks belong to:

for j in $(cat step1_genusspecific_ORFan.txt)
do
    echo "#!/bin/bash" > "$j".sh  # Start each script with a shebang
    echo "for i in \$(grep \"$j\" all_450_proteins.clusters.tsv | cut -f3- -d \"@\" | cut -f1 -d \"(\")" >> "$j".sh
    echo "do" >> "$j".sh
    echo "    echo \"\$i\" > ${j}_temp" >> "$j".sh
    echo "    grep \"\$i\" /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/* | cut -f2- -d \":\" | bedtools sort -i - | bedtools closest -a - -b /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/all_500_gtfs_CDSonly.gtf -D a -id -io | awk -F '\t' '(\$12==\"CDS\")' | cut -f6 -d '\"' >> ${j}_temp" >> "$j".sh
    echo "    grep \"\$i\" /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/* | cut -f2- -d \":\" | bedtools sort -i - | bedtools closest -a - -b /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/all_500_gtfs_CDSonly.gtf -D a -id -io | awk -F '\t' '(\$12==\"CDS\")' | cut -f6 -d '\"' | grep -f - all_450_proteins.clusters.tsv | cut -f1 | sort -u >> ${j}_temp" >> "$j".sh
    echo "    grep \"\$i\" /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/* | cut -f2- -d \":\" | bedtools sort -i - | bedtools closest -a - -b /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/all_500_gtfs_CDSonly.gtf -D a -iu -io | awk -F '\t' '(\$12==\"CDS\")' | cut -f6 -d '\"' >> ${j}_temp" >> "$j".sh
    echo "    grep \"\$i\" /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/* | cut -f2- -d \":\" | bedtools sort -i - | bedtools closest -a - -b /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/all_500_gtfs_CDSonly.gtf -D a -iu -io | awk -F '\t' '(\$12==\"CDS\")' | cut -f6 -d '\"' | grep -f - all_450_proteins.clusters.tsv | cut -f1 | sort -u >> ${j}_temp" >> "$j".sh
    echo "    paste -sd, ${j}_temp | sed \"s/,/\\\t/g\" >> ${j}_flank_analysis" >> "$j".sh
    echo "    grep -E '^([^@]*@){2}[^@]*\$' ${j}_flank_analysis | cut -f3,5 | sort -u > ${j}_flanks" >> "$j".sh
    echo "done" >> "$j".sh
done

ls *.sh | sed "s/^/bash /g" > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh

#Extract gene flanks

for j in $(cat step1_genusspecific_ORFan.txt | cut -f2- -d "@" | cut -f1 -d "(")
do
for i in $(cat -n "$j"_flanks | sed "s/^ *//g" | sed "s/\t/,/g")
do
linenumber=$(echo $i | cut -f1 -d ',')
echo $i | cut -f2 -d ',' | grep --no-group-separator -A1 -f - all_450_CDS.faa | sed "s/>/>"$j"_"$linenumber"_up_/g" >> geneflanks.faa
echo $i | cut -f3 -d ',' | grep --no-group-separator -A1 -f - all_450_CDS.faa | sed "s/>/>"$j"_"$linenumber"_down_/g" >> geneflanks.faa
done
done

#Extract prox flank orthologs:

for j in $(cat step1_genusspecific_ORFan.txt | cut -f2- -d "@" | cut -f1 -d "(")
do
for i in $(grep -E '^([^@]*@){2}[^@]*$' "$j"_flank_analysis | awk -F '\t' '{print $0,$3"_"$5}' | awk '!seen[$6]++' | cat -n | sed "s/^ *//g" | sed "s/\t/,/g")
do
linenumber=$(echo $i | cut -f1 -d ',' )
echo $i | cut -f2 -d ',' | grep -f - /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/all_500_gtfs_CDSonly.gtf | sed "s/ \"/ \""$j"_"$linenumber"_/g" >> proxflanks_interim.gtf
done
done

#Get the sequences
cat proxflanks_interim.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-500,$4,$6,$7,$8,$9}' |  awk -F '\t' '($4>0)' | sed "s/ \"/ \"left500_/g" | gtf2bed | bedtools getfasta -s -name -fi /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/all_500_genomes.fasta -bed - > proxflanks.faa
cat proxflanks_interim.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+500,$6,$7,$8,$9}' | sed "s/ \"/ \"right500_/g" | gtf2bed | bedtools getfasta -s -name -fi /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/all_500_genomes.fasta -bed - >> proxflanks.faa

#Remove those that lack either or both flanks:
grep "^>" proxflanks.faa | sed "s/_/\t/" | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep "," | cut -f 1 -d "(" | grep --no-group-separator -A1 -f - proxflanks.faa > interim
mv interim proxflanks.faa

#Search flanks against genomes:
cd ..
cp flank_tangent/geneflanks.faa .
cp flank_tangent/proxflanks.faa .
makeblastdb -in all_450_genomes.faa -dbtype nucl -out all_450_genomes

blastn -query geneflanks.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out geneflanks_extragenus_blastn
blastn -query geneflanks.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out geneflanks_intragenus_blastn
blastn -query geneflanks.faa -db /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out geneflanks_pangenome_blastn
blastn -query proxflanks.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out proxflanks_extragenus_blastn
blastn -query proxflanks.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out proxflanks_intragenus_blastn
blastn -query proxflanks.faa -db /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out proxflanks_pangenome_blastn

mkdir flanks
cp *_blastn flanks

#targets:

cat proxflanks_extragenus_blastn proxflanks_intragenus_blastn proxflanks_pangenome_blastn | cut -f-2 | sed 's/_/\t/4' | cut -f1,3 | sed "s/_/\t/" | awk '{print$1"\t"$2"%"$3}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | egrep "left.*right|right.*left" | cut -f1 | awk -F'%' '{ values[$1] = (values[$1] == "" ? $2 : values[$1] ", " $2) } END { for (value in values) { print value "\t" values[value] } }' > proxflanks_targets.txt
cat geneflanks_extragenus_blastn geneflanks_intragenus_blastn geneflanks_pangenome_blastn | sed "s/_up_/_up%/g" | sed "s/_down_/_down%/g" | sed "s/\t/%/" | cut -f1 | cut -f1,3 -d "%" | sed "s/%/\t/g" | sed 's/_/\t/3' | awk -F '\t' '{print $2"\t"$1"@"$3}' | sort -u | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | egrep "up.*down|down.*up" | cut -f1 | sed "s/@/%/g" | awk -F'%' '{ values[$1] = (values[$1] == "" ? $2 : values[$1] ", " $2) } END { for (value in values) { print value "\t" values[value] } }' > geneflanks_targets.txt

#targetlist and intervalinfo:

for i in $(cut -f 1 proxflanks_targets.txt)
do
grep -P "$i\t" proxflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 > "$i"_proxflanks_targetlist.txt
cat prox*blastn | grep -w -F -f "$i"_proxflanks_targetlist.txt - | grep "$i"_ | cut -f2- -d "_" | sed 's/_/\t/3' | cut -f1,3- | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > "$i"_proxflanks_intervalinfo
done

#geneflanks:

cat gene*blastn | sort -u > geneflanks_allblastn #for later

for i in $(cut -f 1 geneflanks_targets.txt)
do
grep -P "$i\t" geneflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 > "$i"_geneflanks_targetlist.txt
cat geneflanks_allblastn | grep -w -F -f "$i"_geneflanks_targetlist.txt | grep "$i"_ | sed "s/_up_/\t/g" | sed "s/_down_/\t/g" | cut -f1,3- | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > "$i"_geneflanks_intervalinfo
done

#1. Collapse all varieties of prox and gene flanks into one file per gene cluster:

for i in $(ls *intervalinfo | cut -f1,2 -d "_" | sort -u); do ls "$i"_*intervalinfo | sed "s/^/cat /g" | bash >> "$i"_compiled_intervalinfo.txt; done

#2. Add in names to each intervalinfo file using the file all_contig_protein_taxonomy.tsv, which has been prepared using the code in assigning_conservation_to_genes.sh

for i in $(ls flanks/*_compiled_intervalinfo.txt | cut -f1,2 -d "_" | cut -f2- -d '/' | sort -u)
sort -k1 flanks/"$i"_compiled_intervalinfo.txt -o flanks/"$i"_compiled_intervalinfo.txt
cut -f1 -d " " flanks/"$i"_compiled_intervalinfo.txt | sort -u > temp
grep -w -F -f temp all_contig_protein_taxonomy.tsv | sort -k1 | join -1 1 -2 1 - flanks/"$i"_compiled_intervalinfo.txt | sed 's/ [^ ]*@/ Ecoli@/' > flanks/"$i"_compiled_intervalinfo.taxa.txt
done

#3. Identify presence/absence

#ORFans with flanks:
ls flanks/*_compiled_intervalinfo.txt | cut -f1,2 -d "_" | cut -f2- -d '/' | sort -u > ORFans_w_flanks.txt

#for annotated, intra-genus:
cat all_proteins_vs_noncoliEscherichia_annotated.tsv | awk -F '\t' '($5>60&&$16<0.001)' | grep -F -f ORFans_w_flanks.txt | cut -f1,2 | cut -f 2- -d "@" > intragenus_distribution_interim.txt
#for ORFs, intra-genus:
cat all_proteins_vs_noncoliEscherichia_ORFs_ATG.tsv all_proteins_vs_noncoliEscherichia_ORFs_TTG_GTG.tsv | awk -F '\t' '($5>60&&$16<0.001)' | grep -F -f ORFans_w_flanks.txt | cut -f1,2 | rev | cut -f2- -d "_" | rev | cut -f 2- -d "@" >> intragenus_distribution_interim.txt
sed "s/(+)//g" intragenus_distribution_interim.txt | sed "s/(-)//g" | sort -u | sort -k2 | join -1 2 -2 1 - all_contig_protein_taxonomy.tsv > intragenus_distribution_interim_2.txt

#For pangenome:
cat all_proteins_vs_pangenome_annotated.tsv | awk -F '\t' '($3>80&&$5>70&&$16<0.001)' | cut -f1,2 | sed "s/\t/@/g" | cut -f2,4 | sed "s/@/\t/g" | cut -f2,4 | sed "s/(+)//g" | sed "s/(-)//g" | rev | sed "s/_/@/" | rev | grep -w -F -f ORFans_w_flanks.txt | sed "s/@/_/g" > pangenome_distribution_interim_1.txt
cat all_proteins_vs_pangenome_ORFs.tsv |  awk -F '\t' '($3>80&&$5>70&&$16<0.001)' | cut -f1,2 | cut -f2- -d "@" | grep -w -F -f ORFans_w_flanks.txt | rev | cut -f2- -d "_" | rev | sed "s/(+)//g" | sed "s/(-)//g" | sort -u >> pangenome_distribution_interim_1.txt
sort -u pangenome_distribution_interim_1.txt | sort -k2 | join -1 2 -2 1 - all_contig_protein_taxonomy.tsv > pangenome_distribution_interim_2.txt

#Put them together:
cat intragenus_distribution_interim_2.txt pangenome_distribution_interim_2.txt | cut -f2,3 -d " " | sort -u | sed 's/ [^ ]*@/ Ecoli@/' > presence_absence.interim.tsv

#4. Remove the ones that are present:
#Follow this formula:

for i in $(ls flanks/*_compiled_intervalinfo.taxa.txt | cut -f2- -d '/' | cut -f1,2 -d "_")
do
value=$(grep "$i" /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/all_500_gtfs_CDSonly.gtf | awk -F '\t' '{print $5-$4+1}')
sed -i "s/$/ $value/" flanks/"$i"_compiled_intervalinfo.taxa.txt
grep "$i" presence_absence.interim.tsv | cut -f2 -d " " | grep -v -w -F -f - flanks/"$i"_compiled_intervalinfo.taxa.txt | awk '(($5>($7*0.5))&&($5<10000))' > flanks/"$i"_compiled_intervalinfo.taxa.final.txt
done

#Delete empty files
find . -type f -empty -delete

#To extract sequences, let's compile all individual genomes in the same place:

mv /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoli_Escherichia_individual_genomes/* /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/
mv /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_genomes_individual/* /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/

find /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_genomes_individual/ -type f | xargs -I {} mv {} /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/
find /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoli_Escherichia_individual_genomes/ -type f | xargs -I {} mv {} /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/

#Make DB out of intervals
#DISCOVERED THE MMSEQS

for variable in $(ls flanks/*_compiled_intervalinfo.taxa.final.txt | cut -f2- -d '/' | cut -f1,2 -d "_")
do
    command="grep \"plus\" flanks/${variable}_compiled_intervalinfo.taxa.final.txt | awk '{OFS=\"\"}{print \"samtools faidx /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/\",\$1,\" \",\$1,\":\",\$3,\"-\",\$4}' | sed \"s/\$/ >> flanks\/${variable}_interval.faa/g\" && grep \"minus\" flanks/${variable}_compiled_intervalinfo.taxa.final.txt | awk '{OFS=\"\"}{print \"samtools faidx -i /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/\",\$1,\" \",\$1,\":\",\$3,\"-\",\$4}' | sed \"s/\$/ >> flanks\/${variable}_interval.faa/g\""
    output_file=flanks/"${variable}_joint_samtools.sh"
    bash -c "$command" > "$output_file"
done

ls flanks/*_joint_samtools.sh | sed "s/^/bash /g" > running.sh
conda deactivate
/stor/work/Ochman/hassan/tools/parallelize_run.sh

#Get ORFan CDS sequences for later
grep "^>" step1_genusspecific_ORFans.faa | cut -f2- -d "@" | cut -f1 -d "(" | grep -f - /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/all_500_gtfs_CDSonly.gtf | gtf2bed | bedtools getfasta -s -name -fi /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/all_500_genomes.fasta -bed - > step1_genusspecific_ORFans.CDS.faa

#Still not done:

#Make databases out of intervals
for i in $(ls flanks/*_interval.faa | cut -f2- -d '/' | cut -f1-2 -d "_" | grep -v "HHDFKHIH_03778")
do
makeblastdb -in flanks/"$i"_interval.faa -dbtype nucl -out flanks/"$i"_interval
done


#blastn search
cd flanks
for i in $(ls *_interval.faa | cut -f2- -d '/' | cut -f1-2 -d "_" | grep -v "HHDFKHIH_03778")
do
ls "$i"_interval.faa | cut -f1,2 -d "_" | grep --no-group-separator -A1 -f - /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/step1_genusspecific_ORFans.CDS.faa |
blastn -query - -db "$i"_interval -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out "$i"_interval_blastn -word_size 7
done

#mviewed:

export PERL5LIB=/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/lib/
for i in $(ls flanks/*_interval_blastn | cut -f2- -d '/' | cut -f1,2 -d "_")
do
echo "/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/bin/mview -in blast flanks/${i}_interval_blastn > flanks/${i}_blastn_mviewed"
done

for i in $(ls flanks/*_blastn_mviewed | cut -f2- -d '/' | cut -f1,2 -d "_")
do
querylength=$(grep -A1 "$i" step1_genusspecific_ORFans.CDS.faa | grep -v "^>" | awk '{print length($0)}')
ratio=$(tail -n+8 flanks/"$i"_blastn_mviewed | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}')
if (( $(echo "$ratio > 0.5" | bc -l) ))
then
tail -n+9 flanks/"$i"_blastn_mviewed | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2>70)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear > flanks/"$i"_blastn_seq.faa
tail -n+8 flanks/"$i"_blastn_mviewed | head -n-3 | sed "1s/^/g /g" | awk '{print $2,$NF}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | head -2 >> flanks/"$i"_blastn_seq.faa
fi
done

wc -l *_blastn_seq.faa | grep " 2 " | rev | cut -f 1 -d " " | rev | sed "s/^/rm /g" | bash

#Non-redundancy:
cd flanks
for i in $(ls *_blastn_seq.faa | cut -f1,2 -d "_")
do
grep "^>" "$i"_blastn_seq.faa | cut -f1 -d ":" | tr -d ">" | head -n-1 | sort -u | grep -w -F -f - ../all_contig_protein_taxonomy.tsv | sed 's/\t[^ ]*@/\tEcoli@/' | sort -k1 > interim
seqkit fx2tab "$i"_blastn_seq.faa | sed "s/:/\t/g" | sort -k1 | join -1 1 -2 1 - interim | awk '{print $1,$2, $3":"$4}' | awk '!seen[$3]++' | sed "s/:/ /g" | awk '{print $4":"$1":"$2"\t"$3}' | sed "s/^/>/g" | sed "s/\t/\n/g" > "$i"_mafft_input.faa
tail -2 "$i"_blastn_seq.faa >> "$i"_mafft_input.faa
grep -A1 "$i" ../step1_genusspecific_ORFans.CDS.faa | sed "s/>/>FULL_/g" >> "$i"_mafft_input.faa
done

#mafft step
for i in $(ls *_blastn_seq.faa | cut -f1,2 -d "_")
do
mafft --auto "$i"_mafft_input.faa > "$i"_mafft.aln
done

#To convert to a readable Excel file:

