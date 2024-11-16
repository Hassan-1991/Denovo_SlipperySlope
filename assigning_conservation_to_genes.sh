awk -F '\t' '($12=="Complete Genome")' /stor/scratch/Ochman/hassan/083024_RefSeq_Genbank_genomes/RS_assembly_summary_04062024.txt | cut -f8 | sed "s/\[//g" | sed "s/\]//g" | sed "s/Coxiella-like/Coxiella/g" | grep -v "^unidentified" | grep -v "^endosymbiont" | grep -v "secondary" | sed "s/^uncultured //g" | cut -f1 -d " " | tr -d "\'" | sort -u > genus_names_indatabase.txt

#For annotated proteins, get the taxonomy from the new database
#For the ones missing an assignment in the new database, try older method
grep "^>" /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_proteins/04072024_nonEscherichia_protein.faa | cut -f1 -d " " | tr -d ">" > extragenus_interim1
sort -k1 extragenus_interim1 | join -1 1 -2 2 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/proteins/accession_proteinID_taxonomy.tsv | cut -f1,3 -d " " > extragenus_interim2
#Extract ones that weren't assigned:
cut -f1 -d " " extragenus_interim2 > test
cat extragenus_interim1 test > test2
sort test2 | uniq -c | awk '($1==1)' > test3
#possible shorter method: cut -f1 -d " " extragenus_interim2 | grep -v -F -f - extragenus_interim1 

#ALT:
grep "^>" /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_proteins/04072024_nonEscherichia_protein.faa | cut -f1 -d " " | tr -d ">"

###
grep "^>" /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_proteins/04072024_nonEscherichia_protein.faa > extragenus_interim1 #Extract identifiers from database
sed "s/ /\[/" extragenus_interim1 | awk -F '[' '{print $1,$NF}' | cut -f1-3 -d " " | tr -d ">" | sort -u > extragenus_interim2_annotatedproteins #Massage
cut -f1,2 -d " " extragenus_interim2_annotatedproteins > extragenus_interim2_annotatedproteins_genusonly #Only retain genus information
sort -u extragenus_interim2_annotatedproteins_genusonly -o extragenus_interim2_annotatedproteins_genusonly

#For ORFs:
grep "^>" /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes.getorf.noCTG > extragenus_interim_ORFs_1 #Extract identifiers from database
grep -v "REVERSE SENSE" extragenus_interim_ORFs_1 | cut -f1,5,6 -d " " > extragenus_interim_ORFs_2 #Massage
grep "REVERSE SENSE" extragenus_interim_ORFs_1 | cut -f 1,7,8 -d " " >> extragenus_interim_ORFs_2 #Massage
sed "s/_/ /" extragenus_interim_ORFs_2 | cut -f1,3 -d " " | tr -d ">" | sort -u > extragenus_interim_ORFs_3 #Only retain genus information #code running

#For intragenus (both annotated and ORFs somehow):

grep "seqhdr" /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/gffs/*gff | cut -f 10- -d '/' | sed "s/ /@/" | sed "s/seqhdr=\"/@/g" | cut -f1,3 -d "@" | cut -f1 -d " " | grep "^GCA" | sed 's/_/\t/2' | sed "s/@/\t/g" | awk -F '\t' '{OFS=FS}{print $3,$1}' > Escherichia_contig_accession.tsv
grep "seqhdr" /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/gffs/*gff | cut -f 10- -d '/' | sed "s/ /@/" | sed "s/seqhdr=\"/@/g" | cut -f1,3 -d "@" | cut -f1 -d " " | grep -v "^GCA" | sed "s/\.gff:#@/\t/g" | awk -F '\t' '{OFS=FS}{print $2,$1}' >> Escherichia_contig_accession.tsv
grep ">" /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_proteins.faa  | cut -f1 -d "(" | tr -d ">" | grep -w -F -f - /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/gtfs/*gtf | cut -f1,9 | sed "s/:/\t/g" | cut -f1,3 | sed "s/transcript_id \"//g" | rev | cut -f1 -d "/" | rev | sed "s/transcript_id \"//g" | cut -f1 -d "\"" | cut -f2- -d "/" | sed "s/\.gtf//g" | awk -F '\t' '{OFS=FS}{print $2,$1}' | awk -F'\t' '{if ($2 ~ /GCA_/) sub(/_[^_]+$/, "", $2); print $0}' OFS='\t' | sort -u >> Escherichia_contig_accession.tsv
sed "s/GCA_900636405\.1_41767/GCA_900636405.1/g" Escherichia_contig_accession.tsv | sed "s/GCA_900637015\.1_46514/GCA_900637015\.1/g" | sort -k2 > Escherichia_contig_accession.final.tsv
cat /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/GBRS_all_fastANI.tsv /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/ATB_noncoli_Escherichia_fastANI.tsv | sort -k1  | join -1 1 -2 2 - Escherichia_contig_accession.final.tsv | sed "s/ ATCC 35469//g" | sed "s/ KF1//g" | cut -f3- -d " " | sed "s/Escherichia /Escherichia_/g" | awk '{print $2"\t"$1}' > intragenus_interim_annotated_ORF_1

#For pangenome:
grep "##sequence-region" /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/*gff | grep -v "all_500_gffs" | cut -f 8- -d '/' | sed "s/.gff:##sequence-region /\t/g" | sed "s/ /\t/" | awk -F '\t' '{OFS=FS}{print $2,$1}' > pangenome_contig_interim
grep "^>" all_450_proteins.faa | sed "s/@/\t/g" | tr -d ">"  | cut -f1 -d "(" | awk -F '\t' '{OFS=FS}{print $2,$1}' >> pangenome_contig_interim
sort -u pangenome_contig_interim | sort -k2 > test && mv test pangenome_contig_interim
cut -f1,5 ../../500_ipp_lineagedesignations.tsv | tail -n+2 | sort -k1 | join -1 1 -2 2 - pangenome_contig_interim | sed "s/ /@/" | awk '{print $2"\t"$1}' | sort -k1 > pangenome_contigs_lineages

#Put everything together
cat extragenus_interim2_annotatedproteins_genusonly extragenus_interim_ORFs_3 intragenus_interim_annotated_ORF_1 pangenome_contigs_lineages | sed "s/ /\t/g" > all_contig_protein_taxonomy.tsv
sort -k1 all_contig_protein_taxonomy.tsv -o all_contig_protein_taxonomy.tsv

#Five-step conversion from blast results to distribution levels
#For genera:

#blast to cluster-contig:
cat all_proteins_vs_GBRS_annotated.tsv all_proteins_vs_GBRS_ORFs.tsv | awk -F '\t' '($5>60&&$16<0.001)' | cut -f1,2 | sort -u | awk '{split($2, a, "_"); $2 = a[1]; print}' > cluster_contig.tsv
#cluster-taxa:
sort -k2 cluster_contig.tsv | join -1 2 -2 1 - ../../all_contig_protein_taxonomy.tsv > cluster_contig_taxa.tsv
#cluster-taxa_number:
cut -f2- -d " " cluster_contig_taxa.tsv | sort -u | cut -f1 -d " "  | sort | uniq -c | sed "s/^ *//g" | awk '{print $2,$1}' > cluster_taxanumber.tsv
cut -f1 -d " " cluster_taxanumber.tsv | sed "s/$/\t/g" | grep -vf - ../all_450_proteins.clusters.tsv | cut -f1 | sort -u | sed "s/$/ 0/g" >> cluster_taxanumber.tsv
#gene-cluster-taxa_number:
sort -k1 ../all_450_proteins.clusters.tsv -o ../all_450_proteins.clusters.tsv
sort -k1  cluster_taxanumber.tsv | join -1 1 -2 1 - ../all_450_proteins.clusters.tsv > cluster_taxanumber_gene.tsv
#startposition-gene-cluster-taxa_number:
#For one genome:
awk -F '\t' '($7=="+"&&$3=="CDS")' /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/11657_5#57.gtf | cut -f1,4,9 | cut -f-2 -d "\"" | sed "s/transcript_id \"//g" > temp
awk -F '\t' '($7=="-"&&$3=="CDS")' /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/11657_5#57.gtf | cut -f1,5,9 | cut -f-2 -d "\"" | sed "s/transcript_id \"//g" >> temp
sort -k1 temp -o temp
cut -f1 temp | grep -w -F -f - /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/11657_5#57.gtf | awk -F '\t' '($3=="CDS")' | cut -f1 | sort | uniq -c | sort -nrk1 | sed "s/^ *//g" | sort -k2 > temp2
join -1 1 -2 2 temp temp2 | sed "s/ /\t/g" | sort -k3 > temp3
grep "11657_5#57" cluster_taxanumber_gene.tsv | cut -f2- -d " " | sed "s/ /@/g" | cut -f1 -d "(" | cut -f1,3 -d "@" | sed "s/@/\t/g" | sort -k2 | join -1 2 -2 3 - temp3 | sed "s/ /,/g" | sort -k4,4 -k5,5n | sed "1s/^/gene,frequency,contig,position,order\n/" > cluster_taxanumber_gene_position.tsv

#RETHINK PLOTS
