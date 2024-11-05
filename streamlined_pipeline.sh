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
mmseqs search all_450_proteins all_450_proteins resultDB tmp --min-seq-id 0.8 -c 0.7 --cov-mode 1
mmseqs convertalis all_450_proteins all_450_proteins resultDB resultDB.m8
mmseqs linclust all_450_proteins clusterDB tmp --min-seq-id 0.8 -c 0.7 --cov-mode 1
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
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q clustered_representative_proteins.faa -d all_450_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_annotated.tsv -k 0 -b8 -c1
#Across pangenome, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q clustered_representative_proteins.faa -d all_450_genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_ORFs.tsv -k 0 -b8 -c1

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

#Extract gene flanks #running

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

#
for i in $(ls /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_genomes/*fasta); do seqkit fx2tab $i; done | awk -F '\t' '{OFS=FS}{print $1,length($2)}' >> all_contig_lengths.tsv
sort -k1 all_contig_lengths.tsv -o all_contig_lengths.tsv
cat step1_genusspecific_ORFan.txt | cut -f 2- -d "@" | cut -f1 -d "(" | grep -f - /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/*gtf | cut -f9- -d '/' | sed "s/:/\t/g" | sed "s/\.gtf//g" | cut -f1,2,6,7,11 | cut -f1 -d ";" | tr -d "\"" | sed "s/transcript_id//g" | sort -k2 | join -1 2 -2 1 - all_contig_lengths.tsv > contig_genome_genestart_geneend_genename_contiglength.tsv
cut -f 5 -d " " contig_genome_genestart_geneend_genename_contiglength.tsv | grep -f - all_450_proteins.clusters.tsv | sed "s/\t/@/g" | cut -f1,2,4 -d "@" | rev | sed "s/@/\t/" | rev | sed "s/(+)$//g" | sed "s/(-)$//g" | sort -k2 | join -1 2 -2 5 - contig_genome_genestart_geneend_genename_contiglength.tsv | sed "s/ /\t/g" > genename_cluster_contig_genome_start_end_contiglength.tsv
