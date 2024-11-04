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

for i in $(ls /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_genomes/*fasta); do seqkit fx2tab $i; done | awk -F '\t' '{OFS=FS}{print $1,length($2)}' >> all_contig_lengths.tsv
sort -k1 all_contig_lengths.tsv -o all_contig_lengths.tsv
cat step1_genusspecific_ORFan.txt | cut -f 2- -d "@" | cut -f1 -d "(" | grep -f - /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/*gtf | cut -f9- -d '/' | sed "s/:/\t/g" | sed "s/\.gtf//g" | cut -f1,2,6,7,11 | cut -f1 -d ";" | tr -d "\"" | sed "s/transcript_id//g" | sort -k2 | join -1 2 -2 1 - all_contig_lengths.tsv > contig_genome_genestart_geneend_genename_contiglength.tsv
cut -f 5 -d " " contig_genome_genestart_geneend_genename_contiglength.tsv | grep -f - all_450_proteins.clusters.tsv | sed "s/\t/@/g" | cut -f1,2,4 -d "@" | rev | sed "s/@/\t/" | rev | sed "s/(+)$//g" | sed "s/(-)$//g" | sort -k2 | join -1 2 -2 5 - contig_genome_genestart_geneend_genename_contiglength.tsv | sed "s/ /\t/g" > genename_cluster_contig_genome_start_end_contiglength.tsv



