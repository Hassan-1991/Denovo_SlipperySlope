#Databases for this are stored in ORFan_genomics repository

sed "s/^/>/g" 30547_total | grep --no-group-separator -A1 -F -w -f - F3_pan_genome_reference.prot.fa > 30547_pan_genome_reference.prot.fa

#Search-1: Diamond, vanilla, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q 30547_pan_genome_reference.prot.fa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_proteins/04072024_nonEscherichia_protein.dmnd --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --out Ecoli_30547_vs_GBRS_annotated.tsv -k 0 -b8 -c1
awk -F '\t' '($5>60&&$16<0.001)' Ecoli_30547_vs_GBRS_annotated.tsv | cut -f1 | sort -u > Ecoli_step1_nonORFans.txt
grep -w -F -v -f Ecoli_step1_nonORFans.txt 30547_pan_genome_reference.prot.fa | grep "^>" | tr -d ">" | grep --no-group-separator -w -F -A1 -f - 30547_pan_genome_reference.prot.fa > Ecoli_step1_ORFans.faa

#Search-2: Diamond, ultra-sensitive, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_step1_ORFans.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_proteins/04072024_nonEscherichia_protein.dmnd --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Ecoli_30547_vs_GBRS_annotated_ultrasensitive.tsv -k 0 -b8 -c1
awk -F '\t' '($5>60&&$16<0.001)' Ecoli_30547_vs_GBRS_annotated_ultrasensitive.tsv | cut -f1 | sort -u > Ecoli_step2_nonORFans.txt
grep -w -F -v -f Ecoli_step2_nonORFans.txt Ecoli_step1_ORFans.faa | grep "^>" | tr -d ">" | grep --no-group-separator -w -F -A1 -f - Ecoli_step1_ORFans.faa > Ecoli_step2_ORFans.faa

#Search-3: Diamond, vanilla, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_step2_ORFans.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes.getorf.noCTG.prot.dmnd --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --out Ecolistep2_vs_GBRS_ORFs.tsv -k 0 -b8 -c1
awk -F '\t' '($5>60&&$16<0.001)' Ecolistep2_vs_GBRS_ORFs.tsv | cut -f1 | sort -u > Ecoli_step3_nonORFans.txt
grep -w -F -v -f Ecoli_step3_nonORFans.txt Ecoli_step2_ORFans.faa | grep "^>" | tr -d ">" | grep --no-group-separator -w -F -A1 -f - Ecoli_step2_ORFans.faa > Ecoli_step3_ORFans.faa

#Search-4: Diamond, ultra-sensitive, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_step3_ORFans.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes.getorf.noCTG.prot.dmnd --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Ecolistep2_vs_GBRS_ORFs_ultrasensitive.tsv -k 0 -b8 -c1
awk -F '\t' '($5>60&&$16<0.001)' Ecolistep2_vs_GBRS_ORFs_ultrasensitive.tsv | cut -f1 | sort -u > Ecoli_step4_nonORFans.txt
grep -w -F -v -f Ecoli_step4_nonORFans.txt Ecoli_step3_ORFans.faa | grep "^>" | tr -d ">" | grep --no-group-separator -w -F -A1 -f - Ecoli_step3_ORFans.faa > Ecoli_step4_ORFans.faa

#Flanks:
#Identify the genome (and gene) with the best hit
awk '{if($1 in arr) {if($17 > max[$1]) {arr[$1] = $0; max[$1] = $17}} else {arr[$1] = $0; max[$1] = $17}} END {for(i in arr) print arr[i]}' temp > 30216_best_hit.tsv
grep "^>" Ecoli_step4_ORFans.faa | tr -d ">" | grep -w -F -f - 30216_best_hit.tsv | cut -f1,2 > Ecoli_step4_ORFans.besthits
cut -f 2 -d "@" Ecoli_step4_ORFans.besthits | cut -f1 -d "(" | sed "s/^/=/g" | sed "s/$/;/g" | grep -f - 500_gffs/all_500_gffs.gff > Ecoli_step4_ORFans.besthits.gff
cut -f1 -d ";" Ecoli_step4_ORFans.besthits.gff | sed "s/ID=//g" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' > Ecoli_step4_ORFans.besthits.gtf

cat Ecoli_step4_ORFans.besthits.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-500,$4,$6,$7,$8,$9}' |  awk -F '\t' '($4>0)' | sed "s/ \"/ \"left500_/g" | gtf2bed | bedtools getfasta -s -name -fi 500_gffs/all_500_genomes.fasta -bed - > Ecoli_step4_ORFans_proxflanks.faa
cat Ecoli_step4_ORFans.besthits.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+500,$6,$7,$8,$9}' | sed "s/ \"/ \"right500_/g" | gtf2bed | bedtools getfasta -s -name -fi 500_gffs/all_500_genomes.fasta -bed - >> Ecoli_step4_ORFans_proxflanks.faa

#In cases where the "best hit" genome is on a short contig, I need another candidate
#Get these genes:
grep "^>" Ecoli_step4_ORFans_proxflanks.faa | sed "s/_/\t/" | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep -v "," | cut -f 1 -d "(" > besthits_with_shortcontigs
#270 such cases. How many of them have a second hit?
grep -f besthits_with_shortcontigs Ecoli_step4_ORFans.besthits | cut -f1 | grep -F -w -f - temp | grep -v -F -w -f besthits_with_shortcontigs | awk '{if($1 in arr) {if($17 > max[$1]) {arr[$1] = $0; max[$1] = $17}} else {arr[$1] = $0; max[$1] = $17}} END {for(i in arr) print arr[i]}' | cut -f-2 > 160_shortcontig_besthits
#How many have a second hit on an accommodating contig?
cut -f 2 -d "@" 160_shortcontig_besthits | cut -f1 -d "(" | sed "s/^/=/g" | sed "s/$/;/g" | grep -f - 500_gffs/all_500_gffs.gff > 160shortcontigs.besthits.gff
cut -f1 -d ";" 160shortcontigs.besthits.gff | sed "s/ID=//g" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' > 160shortcontigs.besthits.gtf
cat 160shortcontigs.besthits.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-500,$4,$6,$7,$8,$9}' |  awk -F '\t' '($4>0)' | sed "s/ \"/ \"left500_/g" | gtf2bed | bedtools getfasta -s -name -fi 500_gffs/all_500_genomes.fasta -bed - > 160shortcontigs.besthits.proxflanks.faa
cat 160shortcontigs.besthits.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+500,$6,$7,$8,$9}' | sed "s/ \"/ \"right500_/g" | gtf2bed | bedtools getfasta -s -name -fi 500_gffs/all_500_genomes.fasta -bed - >> 160shortcontigs.besthits.proxflanks.faa
grep "^>" 160shortcontigs.besthits.proxflanks.faa | sed "s/_/\t/" | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep "," | cut -f 1 -d "(" > replacement_besthits
#55 rescued, let's forget about others
#Grounds: the best or second best hit was found on short contig
grep "^>" Ecoli_step4_ORFans_proxflanks.faa | sed "s/_/\t/" | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep "," | cut -f1 -d "(" | grep -f - Ecoli_step4_ORFans.besthits > Ecoli_step4_ORFans.final.besthits
grep -f replacement_besthits 160_shortcontig_besthits >> Ecoli_step4_ORFans.final.besthits
#GTF file
cut -f1 -d "(" Ecoli_step4_ORFans.final.besthits | cut -f 2- -d "@" | grep -f - 500_gffs/500_gtfs_processed/*gtf | cut -f2- -d ":" > Ecoli_step4_ORFans.final.besthits.gtf
#Final proxflanks
grep "^>" Ecoli_step4_ORFans_proxflanks.faa | sed "s/_/\t/" | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep "," | cut -f1 -d "(" | grep -f - Ecoli_step4_ORFans.besthits | cut -f2- -d "@" | grep --no-group-separator -A1 -F -f - Ecoli_step4_ORFans_proxflanks.faa > Ecoli_step4_ORFans.final.proxflanks.faa
grep -f replacement_besthits 160_shortcontig_besthits | cut -f2- -d "@" | grep --no-group-separator -A1 -F -f - 160shortcontigs.besthits.proxflanks.faa >> Ecoli_step4_ORFans.final.proxflanks.faa

#GENEFLANKS
#Make all-encompassing gtf
ls 500_gffs/500_gtfs_processed/ | sed "s/^/cat 500_gffs\/500_gtfs_processed\//g" | bash > 500_gffs/500_gtfs_processed/all_500_gtfs.gtf
for i in $(cut -f2 -d "\"" Ecoli_step4_ORFans.final.besthits.gtf)
do
awk -v var="$i" -F "\"" '($2==var)' Ecoli_step4_ORFans.final.besthits.gtf | bedtools closest -a - -b 500_gffs/500_gtfs_processed/all_500_gtfs.gtf -D a -iu -io | cut -f 10-18 | sed "s/ \"/ \""$i"_down_/g" | sed 's/""/"/g' | sed 's/"_/_/g' | gtf2bed | bedtools getfasta -s -name -fi 500_gffs/all_500_genomes.fasta -bed - >> Ecoli_step4_ORFans_geneflanks.faa
awk -v var="$i" -F "\"" '($2==var)' Ecoli_step4_ORFans.final.besthits.gtf | bedtools closest -a - -b 500_gffs/500_gtfs_processed/all_500_gtfs.gtf -D a -id -io | cut -f 10-18 | sed "s/ \"/ \""$i"_up_/g" | sed 's/""/"/g' | sed 's/"_/_/g' | gtf2bed | bedtools getfasta -s -name -fi 500_gffs/all_500_genomes.fasta -bed - >> Ecoli_step4_ORFans_geneflanks.faa
done

grep "^>" Ecoli_step4_ORFans_geneflanks.faa | cut -f -3 -d "_" | rev | sed "s/_/\t/" | rev | awk '{print $2,$1}' | awk '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep "," | cut -f1 | grep --no-group-separator -A1 -F -f - Ecoli_step4_ORFans_geneflanks.faa > Ecoli_step4_ORFans.final.geneflanks.faa

#blastin time
blastn -query Ecoli_step4_ORFans.final.proxflanks.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Ecoli_step4_ORFans_proxflanks_blastn
cat Ecoli_step4_ORFans_proxflanks_blastn | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | cut -f 1 -d "@" | awk -F '_' '{OFS="\t"}{print $1,$0}' | sed "s/left500_//g" | sed "s/right500_//g" | awk -F '\t' '{OFS=""}{print $1,"\t",$2,"%",$3}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep "left.*right" | cut -f 1 | awk -F'%' '{ values[$1] = (values[$1] == "" ? $2 : values[$1] ", " $2) } END { for (value in values) { print value "\t" values[value] } }' > Ecoli_step4_ORFans_proxflanks_targets.txt
blastn -query Ecoli_step4_ORFans.final.geneflanks.faa -db /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Ecoli_step4_ORFans_geneflanks_blastn
cat Ecoli_step4_ORFans_geneflanks_blastn | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | cut -f 1 -d "@" | sed "s/_up_/%up%/g" | sed "s/_down_/%down%/g" | awk -F "\t" '{OFS=""}{print $1,"%",$2}' | awk -F '%' '{OFS=""}{print $2,"_",$3,"\t",$1,"%",$4}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | egrep "up.*down|down.*up" | cut -f 1 | awk -F'%' '{ values[$1] = (values[$1] == "" ? $2 : values[$1] ", " $2) } END { for (value in values) { print value "\t" values[value] } }' > Ecoli_step4_ORFans_geneflanks_targets.txt

for i in $(cut -f 1 -d "(" Ecoli_step4_ORFans_proxflanks_targets.txt)
do
grep "$i(" Ecoli_step4_ORFans_proxflanks_targets.txt | sed "s/,/\n/g" | sed "s/)\t/)\n/g" | sed "s/^ *//g" | tail -n+2 > "$i"_Ecoli_step4_ORFans_proxflanks_targetlist.txt
grep -f "$i"_Ecoli_step4_ORFans_proxflanks_targetlist.txt Ecoli_step4_ORFans_proxflanks_blastn | grep "$i(" | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | awk -F'\t' 'BEGIN { OFS = "\t" } { sub(/@.*$/, "", $2); print }' | cut -f 2- -d "_" | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > "$i"_Ecoli_step4_ORFans_proxflanks_intervalinfo
done

for i in $(cut -f 1 Ecoli_step4_ORFans_geneflanks_targets.txt)
do
grep -P "$i\t" Ecoli_step4_ORFans_geneflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 > "$i"_Ecoli_step4_ORFans_geneflanks_targetlist.txt
grep -f "$i"_Ecoli_step4_ORFans_geneflanks_targetlist.txt Ecoli_step4_ORFans_geneflanks_blastn | grep "$i"_ | awk -F'\t' -v OFS='\t' '{$2 = $2 "@"}1' | awk -F'\t' 'BEGIN { OFS = "\t" } { sub(/@.*$/, "", $2); print }' | sed "s/_up/%up/g" | sed "s/_down/%down/g" | cut -f 2- -d "%" | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}'  | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' | awk -F '\t' '($4<10000)' > "$i"_Ecoli_step4_ORFans_geneflanks_intervalinfo
done


#What constitutes non-coli Escherichia?

/stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/fastANI_noncoli_Escherichia.sh
#Tag distances with species names
sort -k1 /stor/scratch/Ochman/hassan/0318_AllTheBacteria/hq_dataset.species_calls.tsv | join -1 1 -2 1 - /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/fastANI_noncoli_Esch_interim.tsv | sed "s/Escherichia /Escherichia_/g" | sort -nrk 3 | sed "s/ /\t/g" | awk -F '\t' '{OFS=FS}{print $1,$3,$2}' > /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/ATB_noncoli_Escherichia_fastANI.tsv

#Same with GBRS fastANI
#While in directory /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes
cut -f 2,3 REL606_vs_GBRS_fastANI | sed "s/_genomic.fna//g" | sed "s/_/\t/" | sed "s/_/\t/" | cut -f1,2,4 | sed "s/\t/_/" | sort -k1 > REL606_fastANI_interim
cut -f1,8 /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/GB_assembly_summary_04062024.txt | sort -k1 | join -1 1 -2 1 - REL606_fastANI_interim > GBRS_all_fastANI.tsv
sed "s/ /\t/" GBRS_all_fastANI.tsv | rev | sed "s/ /\t/" | rev | awk -F '\t' '{OFS=FS}{print $1,$3,$2}' | sort -nrk2 > interim
mv interim GBRS_all_fastANI.tsv

cat GBRS_all_fastANI.tsv /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/ATB_noncoli_Escherichia_fastANI.tsv

#After observing the divergences:

#Everything below 93% identity are probably not E. coli
#All reports of E. coli below 93% are problematic, based on manual NCBI searches of the accessions (different warnings)
#This comprises our set of non-coli Escherichia genomes

mkdir 10082024_noncoli_Escherichia_database
cd 10082024_noncoli_Escherichia_database
cat ../GBRS_all_fastANI.tsv /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/ATB_noncoli_Escherichia_fastANI.tsv | sort -nrk 2  | awk -F '\t' '($2<97)' | awk -F '\t' '($2<93)' | awk -F '\t' '($2>84)' | grep -v "Escherichia.*coli" | grep -v "sp" | grep "^GC" | cut -f1 | grep -f - ../../GB_assembly_summary_04062024.txt | cut -f20 | awk -F '/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed "s/^/wget /g" 
cat ../GBRS_all_fastANI.tsv /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/ATB_noncoli_Escherichia_fastANI.tsv | sort -nrk 2  | awk -F '\t' '($2<97)' | awk -F '\t' '($2<93)' | awk -F '\t' '($2>84)' | grep -v "Escherichia.*coli" | grep -v "sp" | grep "^GC" | cut -f1 | grep -f - ../../GB_assembly_summary_04062024.txt | cut -f20 | awk -F '/' '{print $0"/"$NF"_protein.faa.gz"}' | sed "s/^/wget /g" 
cat ../GBRS_all_fastANI.tsv /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/ATB_noncoli_Escherichia_fastANI.tsv | sort -nrk 2  | awk -F '\t' '($2<97)' | awk -F '\t' '($2<93)' | awk -F '\t' '($2>84)' | grep -v "Escherichia.*coli" | grep -v "sp" | grep "^GC" | cut -f1 | grep -f - ../../GB_assembly_summary_04062024.txt | cut -f20 | awk -F '/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed "s/^/wget /g" 

#ATB:
awk -F '\t' '($2<93)' /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/ATB_noncoli_Escherichia_fastANI.tsv | grep -v "Escherichia_sp" | cut -f1 | grep -f - /stor/scratch/Ochman/hassan/0318_AllTheBacteria/sample2species2file.tsv
awk -F '\t' '($2<93)' /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/ATB_noncoli_Escherichia_fastANI.tsv | grep -v "Escherichia_sp" | cut -f1 | grep -f - /stor/scratch/Ochman/hassan/0318_AllTheBacteria/sample2species2file.tsv | awk -F '\t' '{print $3,$1}' | sed "s/.asm.tar.xz /\//g" | sed "s/$/.fa/g" | sed "s/^/\/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\//g" | sed "s/^/cat /g" | bash >> ATB_noncoli_Escherichia.faa
#proteins:
awk -F '\t' '($2<93)' /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/ATB_noncoli_Escherichia_fastANI.tsv | grep -v "Escherichia_sp" | cut -f1 | grep -f - /stor/scratch/Ochman/hassan/0318_AllTheBacteria/sample2species2file.tsv | awk -F '\t' '{print $3,$1}' | sed "s/.asm.tar.xz /\//g" | sed "s/$/.fa/g" | sed "s/^/\/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\//g" | sed "s/^/prodigal -i /g" | sed "s/$/ -f gff -o/g" | sed "s/.fa -f/.fa\/-f/g" | awk -F '/' '{print $0,$9}' | sed "s/.fa\/-f/.fa -f/g" | sed "s/$/t/g" | sed "s/.fat/.gff/g" > running.sh
./parallelize_run.sh

