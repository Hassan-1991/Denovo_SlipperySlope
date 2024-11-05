#For annotated proteins:
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

grep "^" gtfs/*gtf | cut -f1 | rev | cut -f1 -d "/" | rev | sed "s/.gtf:/\t/g" | awk -F '\t' '{OFS=FS}{print $2,$1}' | awk -F'\t' '{if ($2 ~ /GCA_/) sub(/_[^_]+$/, "", $2); print $0}' OFS='\t' | sort -u > Escherichia_contig_accession.tsv
grep ">" /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_proteins.faa  | cut -f1 -d "(" | tr -d ">" | grep -w -F -f - gtfs/*gtf | cut -f1,9 | sed "s/:/\t/g" | cut -f1,3 | sed "s/transcript_id \"//g" | cut -f1 -d "\"" | cut -f2- -d "/" | sed "s/\.gtf//g" | awk -F '\t' '{OFS=FS}{print $2,$1}' | awk -F'\t' '{if ($2 ~ /GCA_/) sub(/_[^_]+$/, "", $2); print $0}' OFS='\t' | sort -u >> Escherichia_contig_accession.tsv
sort -k2 Escherichia_contig_accession.tsv -o Escherichia_contig_accession.tsv
cat /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/GBRS_all_fastANI.tsv /stor/scratch/Ochman/hassan/0318_AllTheBacteria/noncoli_Escherichia/ATB_noncoli_Escherichia_fastANI.tsv | sort -k1  | join -1 1 -2 2 - Escherichia_contig_accession.tsv | sed "s/ ATCC 35469//g" | sed "s/ KF1//g" | cut -f3- -d " " | sed "s/Escherichia /Escherichia_/g" | awk '{print $2"\t"$1}' > intragenus_interim_annotated_ORF_1

#For pangenome:
grep "^" /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/500_gtfs_processed/* | cut -f1 | rev | cut -f1 -d "/" | rev | sed "s/.gtf:/\t/g" | awk -F '\t' '{OFS=FS}{print $2,$1}' > pangenome_contig_interim
grep "^>" all_450_proteins.faa | sed "s/@/\t/g" | tr -d ">"  | cut -f1 -d "(" | awk -F '\t' '{OFS=FS}{print $2,$1}' >> pangenome_contig_interim
sort -u pangenome_contig_interim | sort -k2 > test && mv test pangenome_contig_interim
cut -f1,5 ../../500_ipp_lineagedesignations.tsv | tail -n+2 | sort -k1 | join -1 1 -2 2 - pangenome_contig_interim | sed "s/ /@/" | awk '{print $2"\t"$1}' | sort -k1 > pangenome_contigs_lineages

#Put everything together
cat extragenus_interim2_annotatedproteins_genusonly extragenus_interim_ORFs_3 intragenus_interim_annotated_ORF_1 pangenome_contigs_lineages > all_contig_protein_taxonomy.tsv
sort -k1 all_contig_protein_taxonomy.tsv -o all_contig_protein_taxonomy.tsv
