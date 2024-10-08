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
awk -F '\t' '($5>60&&$16<0.001)' Ecolistep2_vs_GBRS_ORFs.tsv | cut -f1 | sort -u > Ecoli_step3_nonORFans.txt
grep -w -F -v -f Ecoli_step3_nonORFans.txt Ecoli_step2_ORFans.faa | grep "^>" | tr -d ">" | grep --no-group-separator -w -F -A1 -f - Ecoli_step2_ORFans.faa > Ecoli_step3_ORFans.faa



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
