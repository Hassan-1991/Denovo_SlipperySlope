Start with these searches:

#Outside genus, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_proteins/04072024_nonEscherichia_protein --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_GBRS_annotated.tsv -k 0 -b8 -c1
#Outside genus, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes.getorf.noCTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_GBRS_ORFs.tsv -k 0 -b8 -c1
#Outside species, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_annotated.tsv -k 0 -b8 -c1
#Outside species, ORFs_ATG
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes.getorf.ATG --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_ORFs_ATG.tsv -k 0 -b8 -c1
#Outside species, ORFs_TTG_GTG
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes.getorf.TTG_GTG --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_ORFs_TTG_GTG.tsv -k 0 -b8 -c1

#Overall goal:
"""1. Get a taxonomic restriction score (how many species it's present in) for each cluster
2. Assign scores to each gene in the 419 genomes
3. By observing the flanks next to genus-specific ORFans, see if their taxonomic score distribution varies from those of conserved genes
4. See if genome-wide taxonomic restriction scores vary based on within-pangenome lineage
5. Finally, assign each of the genus-specific ORFans to an origin status."""

#Step-1: Identify genus-specific ORFans.
cd /stor/work/Ochman/hassan/Ecoli_pangenome/101224_updated_pipeline
#
grep -w -F -v -f step1_nonORFans.txt ../diamond_searches/all_proteins_reps.faa | grep "^>" | tr -d ">" | grep --no-group-separator -w -F -A1 -f - ../diamond_searches/all_proteins_reps.faa > step1_ORFans.faa
grep -w -F -v -f step1_nonORFans.txt ../diamond_searches/all_proteins_reps.faa | grep "^>" | tr -d ">" > step1_ORFans.txt

awk -F '\t' '($5>50&&$16<0.001)' ../diamond_searches/all_proteins_vs_GBRS_ORFs.tsv > temp
cut -f1 temp | sort -u > ORF_based_genus_nonORFans.txt
grep -Fxvf ORF_based_genus_nonORFans.txt step1_ORFans.txt > step2_ORFans.txt
grep --no-group-separator -w -F -A1 -f step2_ORFans.txt ../diamond_searches/all_proteins_reps.faa > step2_ORFans.faa

#For each cluster, get the genes that are best hits + long contigs

#Get lengths for all contigs
for i in $(ls ../500_gffs/500_genomes/*fasta); do seqkit fx2tab $i; done | awk -F '\t' '{OFS=FS}{print $1,length($2)}' >> contig_lengths.tsv
sort -k1 contig_lengths.tsv -o contig_lengths.tsv
cut -f2 041124_all_protein_clusters.tsv | cut -f 2- -d "@" | cut -f1 -d "(" | grep -f - ../500_gffs/500_gtfs_processed/*gtf | cut -f 4- -d '/' | sed "s/:/\t/g" | sed "s/\.gtf//g" | cut -f1,2,6,7,11 | cut -f1 -d ";" | tr -d "\"" | sed "s/transcript_id//g" | sort -k2 | join -1 2 -2 1 - contig_lengths.tsv > contig_genome_genestart_geneend_genename_contiglength.tsv

#Now identify the "rank" of each protein relative to the representative sequence: ordered by similarity score, which one is the best vs worst match
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d ../500_gffs/419_proteins/all_419_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_reps_vs_419_proteins.tsv -k 0 -b8 -c1
awk -F '\t' '($3>80&&$5>60)' all_reps_vs_419_proteins.tsv | cut -f1,2,17 > cluster_ranking_hits_interim.tsv
awk '{print $0 | "sort -k1,1 -k3,3nr"}' cluster_ranking_hits_interim.tsv | awk 'BEGIN { OFS="\t"; prev_col1=""; rank=0; score=-1 } { if ($1 != prev_col1) {rank=1; score=$3} else if ($3 != score) { rank++; score=$3 } print $0, rank; prev_col1=$1 }' | cut -f1,2,4 | sed "s/\t/%/" | rev | sed "s/@/%/" | rev | cut -f1,3- -d "%" | sed "s/%/\t/g" | rev | sed "s/\t/%/" | rev | sed "s/(+)%/\t/g" | sed "s/(-)%/\t/g" | sort -k2 > cluster_genename_rank.tsv
sort -k5 contig_genome_genestart_geneend_genename_contiglength.tsv | join -1 5 -2 2 - cluster_genename_rank.tsv | sed "s/ /\t/g" > genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv

#Identify best hits

awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==1)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | cut -f7 | sort -u > besthits_1
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==1)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | awk '{print $7,$1}' | awk '!seen[$1]++{print $0}' > all_best_hits.tsv
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==2)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -v -w -f besthits_1 | cut -f7 | sort -u > besthits_2
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==2)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -v -w -f besthits_1 | awk '{print $7,$1}' | awk '!seen[$1]++{print $0}' >> all_best_hits.tsv
cat besthits_1 besthits_2 > besthits_interim
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==3)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -F -v -w -f besthits_interim | cut -f7 | sort -u > temp
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==3)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -F -v -w -f besthits_interim | awk '{print $7,$1}' | awk '!seen[$1]++{print $0}' >> all_best_hits.tsv
cat temp >> besthits_interim
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==4)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -F -v -w -f besthits_interim | cut -f7 | sort -u > temp
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==4)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -F -v -w -f besthits_interim | awk '{print $7,$1}' | awk '!seen[$1]++{print $0}' >> all_best_hits.tsv
cat temp >> besthits_interim
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==5)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -F -v -w -f besthits_interim | cut -f7 | sort -u > temp
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==5)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -F -v -w -f besthits_interim | awk '{print $7,$1}' | awk '!seen[$1]++{print $0}' >> all_best_hits.tsv
cat temp >> besthits_interim
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==6)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -F -v -w -f besthits_interim | cut -f7 | sort -u > temp
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==6)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -F -v -w -f besthits_interim | awk '{print $7,$1}' | awk '!seen[$1]++{print $0}' >> all_best_hits.tsv
cat temp >> besthits_interim
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==7)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -F -v -w -f besthits_interim | cut -f7 | sort -u > temp
awk -F '\t' '((($4-1)>2000)&&(($6-$5)>2000)&&$8==7)' genename_contig_genome_genestart_geneend_contiglength_cluster_rank.tsv | grep -F -v -w -f besthits_interim | awk '{print $7,$1}' | awk '!seen[$1]++{print $0}' >> all_best_hits.tsv
cat temp >> besthits_interim

#Replace representative sequence IDs with actual cluster names
grep "@" besthits_interim | grep -f - 40002_clusters.reps.tsv | cut -f1 | sort -u > besthits_interim_updated
grep -v "@" besthits_interim >> besthits_interim_updated

grep -v -w -f besthits_interim_updated 041124_all_protein_clusters.tsv | cut -f1 | sort -u | grep -w -f - 041124_all_protein_clusters.tsv | cut -f2- -d "@" | cut -f1 -d "(" | grep -w -f - contig_genome_genestart_geneend_genename_contiglength.tsv | awk '((($3-1)>2000)&&(($6-$4)>2000))' | cut -f 5 -d " " | sort -u | grep -w -f - 041124_all_protein_clusters.tsv | awk '!seen[$1]++{print $0}' >> all_best_hits.tsv

#Replace rep seq names with cluster names
sed "s/ /\t/g" all_best_hits.tsv | awk '($1!~"@")' > all_best_hits.final.tsv
sed "s/ /\t/g" all_best_hits.tsv | awk '($1~"@")' | cut -f1 | grep -w -f - 40002_clusters.reps.tsv | cut -f1,2 | sort -u | sort -k2 > replacement
sed "s/ /\t/g" all_best_hits.tsv | awk '($1~"@")' | sort -k1 | join -1 1 -2 2 - replacement  | awk '{print $3"\t"$2}' >> all_best_hits.final.tsv

#Converting step2 ORFans to their best genomic hit:
grep -v "@" step2_ORFans.txt | grep -w -f - all_best_hits.final.tsv > step2_ORFans_besthits.tsv
grep "@" step2_ORFans.txt | grep -f - 40002_clusters.reps.tsv | cut -f1 | sort -u | grep -w -f - all_best_hits.final.tsv >> step2_ORFans_besthits.tsv

#1930(/2314) clusters has at least 1 good sequence per my criteria, after sort -u that comes down to 1860 representative genes.

#get proximal and gene flanks
#First, make a gtf:
cut -f2 step2_ORFans_besthits.tsv | cut -f 1 -d "(" | cut -f2 -d "@" | sort -u | grep -f - ../500_gffs/500_gtfs_processed/*gtf | cut -f2- -d ":" > step2_ORFans_besthits.gtf

#Get 500-bp prox flanks:
cat step2_ORFans_besthits.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-500,$4,$6,$7,$8,$9}' |  awk -F '\t' '($4>0)' | sed "s/ \"/ \"left500_/g" | gtf2bed | bedtools getfasta -s -name -fi ../500_gffs/all_500_genomes.fasta -bed - > step2_ORFans_proxflanks.faa
cat step2_ORFans_besthits.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+500,$6,$7,$8,$9}' | sed "s/ \"/ \"right500_/g" | gtf2bed | bedtools getfasta -s -name -fi ../500_gffs/all_500_genomes.fasta -bed - >> step2_ORFans_proxflanks.faa

#Get gene flanks:
for i in $(cut -f2 -d "\"" step2_ORFans_besthits.gtf)
do
echo "awk -v var=\"${i}\" -F \"\\\"\" '(\$2==var)' step2_ORFans_besthits.gtf | bedtools closest -a - -b ../500_gffs/all_500_gtfs.gtf -D a -iu -io | cut -f 10-18 | sed \"s/ \\\"/ \\\""${i}"_down_/g\" | sed \"s/\\\"\\\"/\\\"/g\" | sed \"s/\\\"_/_/g\" | gtf2bed | bedtools getfasta -s -name -fi ../500_gffs/all_500_genomes.fasta -bed - >> step2_ORFans_geneflanks.faa"
echo "awk -v var=\"${i}\" -F \"\\\"\" '(\$2==var)' step2_ORFans_besthits.gtf | bedtools closest -a - -b ../500_gffs/all_500_gtfs.gtf -D a -id -io | cut -f 10-18 | sed \"s/ \\\"/ \\\""${i}"_up_/g\" | sed \"s/\\\"\\\"/\\\"/g\" | sed \"s/\\\"_/_/g\" | gtf2bed | bedtools getfasta -s -name -fi ../500_gffs/all_500_genomes.fasta -bed - >> step2_ORFans_geneflanks.faa"
done | awk '{print $0"_"int((NR-1)/40)+1}' | split -l 40

ls x* | sed "s/^/bash /g" > running.sh
./parallelize_run.sh

cat step2_ORFans_geneflanks.faa_* > step2_ORFans_geneflanks.faa

#Prepare the best hit sequences (DNA and protein) for later steps:
cat step2_ORFans_besthits.gtf | gtf2bed | bedtools getfasta -s -name -fi ../500_gffs/all_500_genomes.fasta -bed - > step2_ORFans_wflanks_besthits.faa
/stor/work/Ochman/hassan/tools/faTrans step2_ORFans_wflanks_besthits.faa step2_ORFans_wflanks_besthits.prot.faa

#Flank part begins

blastn -query step2_ORFans_proxflanks.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out step2_ORFans_proxflanks_blastn
cat step2_ORFans_proxflanks_blastn | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | cut -f 1 -d "@" | awk -F '_' '{OFS="\t"}{print $1,$0}' | sed "s/left500_//g" | sed "s/right500_//g" | awk -F '\t' '{OFS=""}{print $1,"\t",$2,"%",$3}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep "left.*right" | cut -f 1 | awk -F'%' '{ values[$1] = (values[$1] == "" ? $2 : values[$1] ", " $2) } END { for (value in values) { print value "\t" values[value] } }' > step2_ORFans_proxflanks_targets.txt
blastn -query step2_ORFans_geneflanks.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out step2_ORFans_geneflanks_blastn
cat step2_ORFans_geneflanks_blastn | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | cut -f 1 -d "@" | sed "s/_up_/%up%/g" | sed "s/_down_/%down%/g" | awk -F "\t" '{OFS=""}{print $1,"%",$2}' | awk -F '%' '{OFS=""}{print $2,"_",$3,"\t",$1,"%",$4}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | egrep "up.*down|down.*up" | cut -f 1 | awk -F'%' '{ values[$1] = (values[$1] == "" ? $2 : values[$1] ", " $2) } END { for (value in values) { print value "\t" values[value] } }' > step2_ORFans_geneflanks_targets.txt

for i in $(cut -f 1 -d "(" step2_ORFans_proxflanks_targets.txt)
do
grep "$i(" step2_ORFans_proxflanks_targets.txt | sed "s/,/\n/g" | sed "s/)\t/)\n/g" | sed "s/^ *//g" | tail -n+2 > "$i"_step2_ORFans_proxflanks_targetlist.txt
grep -f "$i"_step2_ORFans_proxflanks_targetlist.txt step2_ORFans_proxflanks_blastn | grep "$i(" | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | awk -F'\t' 'BEGIN { OFS = "\t" } { sub(/@.*$/, "", $2); print }' | cut -f 2- -d "_" | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > "$i"_step2_ORFans_proxflanks_intervalinfo
done

for i in $(cut -f 1 step2_ORFans_geneflanks_targets.txt)
do
grep -P "$i\t" step2_ORFans_geneflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 > "$i"_step2_ORFans_geneflanks_targetlist.txt
grep -f "$i"_step2_ORFans_geneflanks_targetlist.txt step2_ORFans_geneflanks_blastn | grep "$i"_ | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | awk -F'\t' 'BEGIN { OFS = "\t" } { sub(/@.*$/, "", $2); print }' | sed "s/_up/%up/g" | sed "s/_down/%down/g" | cut -f 2- -d "%" | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}'  | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' | awk -F '\t' '($4<10000)' > "$i"_step2_ORFans_geneflanks_intervalinfo
done

for i in $(ls *prox*info | cut -f1,2 -d "_")
do
    value2=$(grep "$i" step2_ORFans_besthits.gtf | awk -F '\t' '{print $5-$4+1}')
    sed -i "s/$/ $value2/" "$i"_step2_ORFans_proxflanks_intervalinfo
done

for i in $(ls *gene*info | cut -f1,2 -d "_")
do
    value2=$(grep "$i" step2_ORFans_besthits.gtf | awk -F '\t' '{print $5-$4+1}')
    sed -i "s/$/ $value2/" "$i"_step2_ORFans_geneflanks_intervalinfo
done

#sidebar: meaningful flanks for comparative analysis between ORFans and non-ORFans
for i in $(ls *prox*info | cut -f1,2 -d "_"); do awk '(($4<($6+600))&&($4>($6*0.5)))' "$i"_step2_ORFans_proxflanks_intervalinfo > "$i"_step2_ORFans_proxflanks_intervalinfo_2; done
for i in $(ls *gene*info | cut -f1,2 -d "_"); do awk '(($4<($6+2000))&&($4>($6*0.5)))' "$i"_step2_ORFans_geneflanks_intervalinfo > "$i"_step2_ORFans_geneflanks_intervalinfo_2; done
find . -type f -empty -delete

#At this point, the two streams can hopefully join

for i in $(ls *intervalinfo_2 | cut -f1,2 -d "_" | sort -u); do ls *"$i"_*intervalinfo_2 | sed "s/^/cat /g" | bash >> "$i"_joint_intervalcoords; done

#Make DB out of intervals

for variable in $(ls *_joint_intervalcoords | rev | cut -f 3- -d "_" | rev); do
    command="grep \"plus\" ${variable}_joint_intervalcoords | awk '{OFS=\"\"}{print \"samtools faidx /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/\",\$1,\" \",\$1,\":\",\$2,\"-\",\$3}' | sed \"s/\$/ >> ${variable}_interval.faa/g\" && grep \"minus\" ${variable}_joint_intervalcoords | awk '{OFS=\"\"}{print \"samtools faidx -i /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/individual_genomes/\",\$1,\" \",\$1,\":\",\$2,\"-\",\$3}' | sed \"s/\$/ >> ${variable}_interval.faa/g\""
    output_file="${variable}_joint_samtools.sh"
    bash -c "$command" > "$output_file"
done

ls *_joint_samtools.sh | sed "s/^/bash /g" > running_sam.sh

#Make databases out of intervals
for i in $(ls *_interval.faa | cut -f1,2 -d "_")
do
makeblastdb -in "$i"_interval.faa -dbtype nucl -out "$i"_interval
done

#tblastn
for i in $(ls *_interval.faa | cut -f1,2 -d "_")
do
grep --no-group-separator -A1 "$i" ../step2_ORFans_wflanks_besthits.prot.faa |
tblastn -query - -db "$i"_interval -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out "$i"_interval_tblastn
done

#blastn
for i in $(ls *_interval.faa | cut -f1,2 -d "_")
do
grep --no-group-separator -A1 "$i(" ../step2_ORFans_wflanks_besthits.faa |
blastn -query - -db "$i"_interval -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out "$i"_interval_blastn -word_size 7
done

export PERL5LIB=/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/lib/

#tblastn, contd
for i in $(ls *_interval_tblastn | cut -f1,2 -d "_")
do
    echo "/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/bin/mview -in blast ${i}_interval_tblastn > ${i}_tblastn_mviewed"
done
for i in $(ls *tblastn_mviewed | cut -f1,2 -d "_")
do
querylength=$(grep -A1 "$i(" ../step2_ORFans_wflanks_besthits.prot.faa | grep -v "^>" | awk '{print length($0)}')
ratio=$(tail -n+7 "$i"_tblastn_mviewed | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}')
if (( $(echo "$ratio > 0.5" | bc -l) ))
then
tail -n+8 "$i"_tblastn_mviewed | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2>70)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | seqkit rmdup -s | linear > "$i"_tblastn_seq.aln
tail -n+7 "$i"_tblastn_mviewed | head -n-3 | sed "1s/^/g /g" | awk '{print $2,$NF}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | head -2 >> "$i"_tblastn_seq.aln
grep -A1 "$i(" ../step2_ORFans_wflanks_besthits.prot.faa | sed "s/>/>FULL_/g" >> "$i"_tblastn_seq.aln
fi
done
#Remove ones that don't meet this criteria
wc -l *_tblastn_seq.aln | grep " 4 " | rev | cut -f 1 -d " " | rev | sed "s/^/rm /g" | bash

#blastn, contd
for i in $(ls *_interval_blastn | cut -f1,2 -d "_")
do
echo "/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/bin/mview -in blast ${i}_interval_blastn > ${i}_blastn_mviewed"
done
for i in $(ls *_blastn_mviewed | cut -f1,2 -d "_")
do
querylength=$(grep -A1 "$i(" ../step2_ORFans_wflanks_besthits.faa | grep -v "^>" | awk '{print length($0)}')
ratio=$(tail -n+8 "$i"_blastn_mviewed | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}')
if (( $(echo "$ratio > 0.5" | bc -l) ))
then
tail -n+9 "$i"_blastn_mviewed | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2>70)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | seqkit rmdup -s | linear > "$i"_blastn_seq.faa
tail -n+8 "$i"_blastn_mviewed | head -n-3 | sed "1s/^/g /g" | awk '{print $2,$NF}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | head -2 >> "$i"_blastn_seq.faa
grep -A1 "$i(" ../step2_ORFans_wflanks_besthits.faa | sed "s/>/>FULL_/g" >> "$i"_blastn_seq.faa
fi
done
wc -l *_blastn_seq.faa | grep " 4 " | rev | cut -f 1 -d " " | rev | sed "s/^/rm /g" | bash

#macse
ls *_blastn_seq.faa | awk '{OFS=""}{print "java -jar \/stor\/work\/Ochman\/hassan\/tools\/macse_v2.06.jar -prog alignSequences -fs 1000 -fs_term 1000 -stop 100 -seq ",$1," -out_NT ",$1,"_NT.aln -out_AA ",$1,"_AA.aln"}' > running.sh

#For tallying purposes, I should make a tsv with the cluster and their best hits for genus-specific ORFans
cut -f2 -d "\"" step2_ORFans_besthits.gtf | grep -w -f - all_best_hits.final.tsv | sed "s/\t/%/g" | sed "s/@/%/g" | egrep ")$" | cut -f1,3 -d "%" | sed "s/%/\t/g" | cut -f1 -d "(" > step2_ORFans_besthits.tsv
cut -f2 -d "\"" step2_ORFans_besthits.gtf | grep -w -f - all_best_hits.final.tsv | grep -v ")$" >> step2_ORFans_besthits.tsv
awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' step2_ORFans_besthits.tsv > step2_ORFans_besthits.final.tsv


#########

#Intra-genus searches
#Fill out an excel sheet with present non-coding, no flanks and flanks states

cat all_proteins_vs_noncoliEscherichia_annotated.tsv all_proteins_vs_noncoliEscherichia_ORFs_TTG_GTG.tsv all_proteins_vs_noncoliEscherichia_ORFs_ATG.tsv | grep -w -f step2_ORFans_besthits.repseqversion.txt > step2_ORFans_besthits.noncoliEscherichia_hits_interim.tsv
awk -F '\t' '($5>60&&$16<0.001)' step2_ORFans_besthits.noncoliEscherichia_hits_interim.tsv | cut -f1,2 | sort -u > step2_ORFans_besthits.noncoliEscherichia_hits.tsv
sed "s/ /_/g" noncoliEscherichia_genomes_ANI_species_contigs.tsv | sort -k4 | join -1 4 -2 2 - step2_ORFans_besthits.noncoliEscherichia_hits.tsv_interim | awk '{OFS=FS}{print $NF,$(NF-1)}' | sed "s/Escherichia /Escherichia_/g" | sed "s/ /\t/" | sed "s/_KF1//g" | sed "s/_ATCC_35469//g" | sort -u | awk '{print $2,$1}' | sed "s/ /\t/g" | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' > Escherichia_presence_absence_interim.tsv

awk '{
    fergusonii = ($0 ~ /fergusonii/) ? 1 : 0;
    ruysiae = ($0 ~ /ruysiae/) ? 1 : 0;
    marmotae = ($0 ~ /marmotae/) ? 1 : 0;
    whittamii = ($0 ~ /whittamii/) ? 1 : 0;
    albertii = ($0 ~ /albertii/) ? 1 : 0;
    print $1, fergusonii "," ruysiae "," marmotae "," whittamii "," albertii
}' Escherichia_presence_absence_interim.tsv > Escherichia_presence_absence.tsv

cut -f1 -d " " Escherichia_presence_absence.tsv | grep -v -w -f - step2_ORFans_besthits.repseqversion.txt | sed "s/$/ 0,0,0,0,0/g" >> Escherichia_presence_absence.tsv

#Next steps: think about how to fill in the 0s
#Will have to tally meaningful flank presence/absence similarly as well

#SIDEBAR
#For assigning taxonomic restriction to all genes:
cat ../diamond_searches/all_proteins_vs_GBRS_annotated.tsv_xaa ../diamond_searches/all_proteins_vs_GBRS_annotated.tsv_xab ../diamond_searches/all_proteins_vs_GBRS_annotated.tsv_xac ../diamond_searches/all_proteins_vs_GBRS_annotated.tsv_xad | awk -F '\t' '($5>60&&$16<0.001)' | cut -f1,2 > all_nonEscherichia_hits_query_targets.tsv
cat ../diamond_searches/all_proteins_vs_noncoliEscherichia_annotated.tsv ../diamond_searches/all_proteins_vs_noncoliEscherichia_ORFs_ATG.tsv ../diamond_searches/all_proteins_vs_noncoliEscherichia_ORFs_TTG_GTG.tsv | awk -F '\t' '($5>60&&$16<0.001)' | cut -f1,2 > all_noncoli_Escherichia_hits_query_targets.tsv

#Getting taxonomic information for each hit:
cd /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_proteins/
#For proteins:
for i in $(ls | grep "GC.*faa" | cut -f1,2 -d "_")
do
grep "^>" "$i"_*_protein.faa | sed "s/^/"$i"\t/g" >> proteins_filenames_contignames_interim.tsv
done
#For genomes and ORFs:
grep "^>" /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes.faa > genomes_filenames_contignames_interim.tsv
cut -f1,2,3 -d " " genomes_filenames_contignames_interim.tsv | tr -d ">" | sort -u | sed "s/ /\t/" > nonEscherichia_genomes_contigs_taxonomy.tsv
sort -k1 nonEscherichia_genomes_contigs_taxonomy.tsv -o nonEscherichia_genomes_contigs_taxonomy.tsv
#Hope this works...let's see if these cover all the "hit" IDs
#I need the taxonomy associated with each contig ID
#For genus outgroups:
