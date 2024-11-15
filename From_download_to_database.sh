#This code shows how to construct the genome and protein databases required to detect de novo genes
#Using E. coli, Salmonella enterica serovar Typhimurium, and Mycobacterium tuberculosis as examples
#Pangenome databases are not included, since there's still no consistent way of making a species pangenome
#But would possibly be added later
#Also, the ATB database download step is also not shown, although that's a simple issue

#Download GenBank assembly summary:

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
mv assembly_summary.txt 102824_GenBank_bacteria_assembly_summary.txt
awk -F '\t' '($12=="Complete Genome")' 102824_GenBank_bacteria_assembly_summary.txt | cut -f20 | awk -F '/' '{print $0"/"$NF"_genomic.fna.gz"}' > 102824_GenBank_bacteria_genome_URLs.txt
awk -F '\t' '($12=="Complete Genome")' 102824_GenBank_bacteria_assembly_summary.txt | cut -f20 | awk -F '/' '{print $0"/"$NF"_protein.faa.gz"}' > 102824_GenBank_bacteria_protein_URLs.txt
mkdir genomes && cp 102824_GenBank_bacteria_genome_URLs.txt genomes
cd genomes && wget -i 102824_GenBank_bacteria_genome_URLs.txt #do this in screen
mkdir proteins && cp 102824_GenBank_bacteria_protein_URLs.txt proteins
cd proteins && wget -i 102824_GenBank_bacteria_protein_URLs.txt

#After downloading and gunzipping, make a table of accession/file name - species - contig name
#For accessions and taxonomy:
awk -F '\t' '($12=="Complete Genome")' 102824_GenBank_bacteria_assembly_summary.txt | cut -f1,8 | sed "s/ /_/g" > accession_taxonomy.tsv #retaining full taxonomy ID for now
sort -k1 accession_taxonomy.tsv -o accession_taxonomy.tsv
#For proteins and accessions:
cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes/proteins
for i in $(ls | grep "GCA" | cut -f 1,2 -d "_"); do grep "^>" "$i"*protein.faa | cut -f1 -d " " | tr -d ">" | sed "s/^/"$i"\t/g" >> accession_proteinID.tsv; done
#tie this to taxonomy:
sort -k1 accession_proteinID.tsv | join -1 1 -2 1 - ../accession_taxonomy.tsv | sed "s/ /\t/g" > accession_proteinID_taxonomy.tsv
sort -k2 /stor/scratch/Ochman/hassan/100724_Complete_Genomes/proteins/accession_proteinID_taxonomy.tsv -o /stor/scratch/Ochman/hassan/100724_Complete_Genomes/proteins/accession_proteinID_taxonomy.tsv

#Make ORF database from all genomes:
for i in $(ls | grep "GCA" | cut -f 1,2 -d "_"); do grep "^>" "$i"*genomic.fna | cut -f1 -d " " | tr -d ">" | sed "s/^/"$i"\t/g" >> accession_genomeID.tsv; done
sort -k1 accession_genomeID.tsv | join -1 1 -2 1 - ../accession_taxonomy.tsv | sed "s/ /\t/g" > accession_genomeID_taxonomy.tsv

#Finally, all the genomes need to dumped in one directory for access at a later step:
#Making the individual genomes folder:
#Concatenate all GenBank genomes, all species external genus genomes in one file, linearize, split -l2, rename
#For ATB:
egrep "Escherichia|Salmonella|Mycobacterium" /stor/scratch/Ochman/hassan/0318_AllTheBacteria/hq_dataset.species_calls.tsv | egrep -v "Escherichia coli|Mycobacterium tuberculosis|Salmonella enterica" | cut -f1 | sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\/\*\//g" | sed "s/$/.fa >> all_ATB_genomes.faa/g" | bash
#For the rest:
ls /stor/scratch/Ochman/hassan/100724_Complete_Genomes/genomes/ | grep "GC.*fna" | sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" | sed "s/$/ >> all_GB_genomes.faa/g" | bash
#Once done:
cat all_ATB_genomes.faa all_GB_genomes.faa | seqkit fx2tab | sed "s/^/>/g" | sed "s/\t/\n/g" | sed '/^$/d' | split -l 2 #each x-file contains one contig and corresponding sequence
for i in x* #renaming by the name of contig
do
echo "mv ${i} \"\$(grep \"^>\" ${i} | cut -f1 -d \" \" | tr -d \">\")\""
done

#To make sure these genomes are what it says on the tin, we calculate ANI from one focal strain in each species to all downloaded genomes

#E. coli focal strain: /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/sequence_RS.fasta
#S. enterica focal strain: /stor/work/Ochman/hassan/protogene_extension/Salmonella_list/GCA_000210855.2_ASM21085v2_genomic.fna
#M. tuberculosis focal strain: /stor/work/Ochman/hassan/protogene_extension/Mycobacterium_list/H37Rv.fna

#######################################################################################################################################
#Escherichia:
cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes/
mkdir Escherichia_temp && cd Escherichia_temp
cp /stor/work/Ochman/hassan/tools/parallelize_run.sh .
#For GB:
ls /stor/scratch/Ochman/hassan/100724_Complete_Genomes/genomes | grep "^GC" | cut -f1,2 -d "_" | sed "s/^/fastANI -q \/stor\/work\/Ochman\/hassan\/protogene_extension\/Ecoli_list\/sequence_RS.fasta -r \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" | sed "s/$/\*fna/g" | awk -F '/' '{print $0" -o "$NF}' | sed "s/$/t/g" | sed "s/\*fnat/_fastANI/g" | split -l 600
ls x* | sed "s/^/bash /g" > running.sh
./parallelize_run.sh
#For ATB: (could be made more non-redundant)
rm x*
egrep "Escherichia" /stor/scratch/Ochman/hassan/0318_AllTheBacteria/hq_dataset.species_calls.tsv | grep -v "Escherichia coli" | cut -f1 | sed "s/^/fastANI -q \/stor\/work\/Ochman\/hassan\/protogene_extension\/Ecoli_list\/sequence_RS.fasta -r \/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\/\*\//g" | sed "s/$/\.fa/g" | awk -F '/' '{print $0" -o "$NF}' | sed "s/$/t/g" | sed "s/\.fat/_fastANI/g" | split -l 100
ls x* | sed "s/^/bash /g" > running.sh
./parallelize_run.sh

ls | egrep "fastANI$" | sed "s/^/cat /g" | sed "s/$/ >> Escherichia_all_fastANI/g" | bash #To prevent argument list too long error

awk -F '\t' '($2~"GCA_")' Escherichia_all_fastANI | cut -f2,3 | sed "s/\/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\///g" | sed 's/_/\t/2' | cut -f1,3 > Escherichia_accession_ANI.tsv
awk -F '\t' '($2!~"GCA_")' Escherichia_all_fastANI | cut -f2,3 | rev | cut -f1 -d '/' | rev | sed "s/\.fa\t/\t/g" >> Escherichia_accession_ANI.tsv
sort -k1 Escherichia_accession_ANI.tsv -o Escherichia_accession_ANI.tsv

#Add taxonomy:
grep "^GCA" Escherichia_accession_ANI.tsv | cut -f1 | grep -w -F -f - ../102824_GenBank_bacteria_assembly_summary.txt | sort -u | cut -f1,8 | sed "s/ /_/g" | sort -k1 | join -1 1 -2 1 - Escherichia_accession_ANI.tsv > Escherichia_taxonomy_accession_ANI.tsv
grep -v "^GCA" Escherichia_accession_ANI.tsv | cut -f1 | grep -w -F -f - /stor/scratch/Ochman/hassan/0318_AllTheBacteria/hq_dataset.species_calls.tsv | sed "s/ /_/g" | sort -k1 | join -1 1 -2 1 - Escherichia_accession_ANI.tsv >> Escherichia_taxonomy_accession_ANI.tsv
sort -nrk3 Escherichia_taxonomy_accession_ANI.tsv | sed "s/ /\t/g" > test && mv test Escherichia_taxonomy_accession_ANI.tsv

#Genus-excluded accessions:
awk -F '\t' '($3>89)' Escherichia_taxonomy_accession_ANI.tsv | grep "GCA_" | cut -f1 > Escherichia_genusexcluded_accesions.txt
awk -F '\t' '($3<89)' Escherichia_taxonomy_accession_ANI.tsv | egrep "Escherichia|Shigella" | cut -f1 >> Escherichia_genusexcluded_accesions.txt

cd ../Escherichia_db
#Escherichia_excluded genomes:
ls ../genomes/ | grep "^GCA" | grep -F -v -f ../Escherichia_temp/Escherichia_genusexcluded_accesions.txt | sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" | sed "s/$/ >> Escherichia_excluded.genomes.faa/g" | bash
#Escherichia_excluded proteins:
ls ../proteins/ | grep "^GCA" | grep -F -v -f ../Escherichia_temp/Escherichia_genusexcluded_accesions.txt | sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/proteins\//g" | sed "s/$/ >> Escherichia_excluded.proteins.faa/g" | bash

#Escherichia_excluded_ORFs:
getorf -sequence Escherichia_excluded.genomes.faa -outseq Escherichia_excluded.genomes.getorf.all -table 1 -minsize 30 -find 3
seqkit fx2tab Escherichia_excluded.genomes.getorf.all | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/" > Escherichia_excluded.genomes.getorf.ATG_TTG_GTG
/stor/work/Ochman/hassan/tools/faTrans -stop Escherichia_excluded.genomes.getorf.ATG_TTG_GTG Escherichia_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa
pwd

#Ecoli_excluded genomes:
awk -F '\t' '($3<93&&$3>89)' ../Escherichia_temp/Escherichia_taxonomy_accession_ANI.tsv | egrep -v "Escherichia_coli|Escherichia_sp" | grep "^GCA" | cut -f1 > Ecoli_excluded_accessions.GBRS.txt
awk -F '\t' '($3<93&&$3>89)' ../Escherichia_temp/Escherichia_taxonomy_accession_ANI.tsv | egrep -v "Escherichia_coli|Escherichia_sp" | grep -v "^GCA" | cut -f1 > Ecoli_excluded_accessions.ATB.txt
ls ../genomes/ | grep "^GCA" | grep -F -f Ecoli_excluded_accessions.GBRS.txt | sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" | sed "s/$/ >> Ecoli_excluded.genomes.faa/g" | bash
sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" Ecoli_excluded_accessions.GBRS.txt | sed "s/$/* >> Ecoli_excluded.genomes.faa/g" | bash
sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\/\*\//g" Ecoli_excluded_accessions.ATB.txt | sed "s/$/\.fa >> Ecoli_excluded.genomes.faa/g" | bash

#Ecoli_excluded_proteins:
sed "s/^/ls \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" Ecoli_excluded_accessions.GBRS.txt | sed "s/$/*/g" | bash > prodigal_interim
sed "s/^/ls \/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\/\*\//g" Ecoli_excluded_accessions.ATB.txt | sed "s/$/\.fa/g" | bash >> prodigal_interim

#Run prodigal:
sed "s/^/prodigal -i /g" prodigal_interim | sed "s/$/ -f gff -o /g" | awk -F '/' '{print $0$NF}' | sed "s/ -f gff -o $//g" | sed "s/_genomic\.fna$/\.gff/g" | sed "s/\.fa$/\.gff/g" | split -l100
ls x* | sed "s/^/bash /g" > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh
#gff to gtf:
for i in $(ls *gff | rev | cut -f2- -d "." | rev); do grep -v "#" "$i".gff | awk -F '\t' '($3=="CDS")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,$1}' | rev | cut -f 2- -d "." | rev | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,"transcript_id \""$9"_"NR"\";gene_id \""$9"_"NR"\";"}'; done

#make space for the newly generated files:
mkdir gff
mkdir gtf
mkdir cds
mkdir proteins
mv *.gff gff
mv *.gtf gtf

#gtf to cds:
for i in $(ls gtf/GCA*gtf | rev | cut -f2- -d "." | rev | cut -f2- -d "/"); do cat gtf/"$i".gtf | gtf2bed | bedtools getfasta -s -name -fi /stor/scratch/Ochman/hassan/100724_Complete_Genomes/genomes/"$i"_genomic.fna -bed - > cds/"$i".cds.faa; done
for i in $(ls gtf/SA*gtf | rev | cut -f2- -d "." | rev | cut -f2- -d "/"); do cat gtf/"$i".gtf | gtf2bed | bedtools getfasta -s -name -fi /stor/scratch/Ochman/hassan/0318_AllTheBacteria/AllTheBacteria_OG/*/$i.fa -bed - > cds/$i.cds.faa; done

#cds to protein:
for i in $(ls cds/*faa | rev | cut -f3- -d '.' | rev | cut -f2- -d "/"); do /stor/work/Ochman/hassan/tools/faTrans -stop cds/$i.cds.faa proteins/$i.protein.faa; done

#Ecoli_excluded_proteins:
cat proteins/*faa >> Ecoli_excluded.proteins.faa

#Ecoli_excluded_ORFs:
getorf -sequence Ecoli_excluded.genomes.faa -outseq Ecoli_excluded.genomes.getorf.all -table 1 -minsize 30 -find 3
seqkit fx2tab Ecoli_excluded.genomes.getorf.all | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/" > Ecoli_excluded.genomes.getorf.ATG_TTG_GTG
/stor/work/Ochman/hassan/tools/faTrans -stop Ecoli_excluded.genomes.getorf.ATG_TTG_GTG Ecoli_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa

#Make diamond and blastn databases:

makeblastdb -in Escherichia_excluded.genomes.faa -dbtype nucl -out Escherichia_excluded.genomes
makeblastdb -in Ecoli_excluded.genomes.faa -dbtype nucl -out Ecoli_excluded.genomes
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Escherichia_excluded.proteins.faa --db Escherichia_excluded.proteins
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Ecoli_excluded.proteins.faa --db Ecoli_excluded.proteins
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Escherichia_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa --db Escherichia_excluded.genomes.getorf.ATG_TTG_GTG.prot
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Ecoli_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa --db Ecoli_excluded.genomes.getorf.ATG_TTG_GTG.prot

#######################################################################################################################################
#Salmonella:
cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes/
mkdir Salmonella_temp && cd Salmonella_temp
cp /stor/work/Ochman/hassan/tools/parallelize_run.sh .
ls /stor/scratch/Ochman/hassan/100724_Complete_Genomes/genomes | grep "^GC" | cut -f1,2 -d "_" | sed "s/^/fastANI -q \/stor\/work\/Ochman\/hassan\/protogene_extension\/Salmonella_list\/GCA_000210855.2_ASM21085v2_genomic.fna -r \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" | sed "s/$/\*fna/g" | awk -F '/' '{print $0" -o "$NF}' | sed "s/$/t/g" | sed "s/\*fnat/_fastANI/g" | split -l 600
ls x* | sed "s/^/bash /g" > running.sh
./parallelize_run.sh
#For ATB: (could be made more non-redundant)
rm x*
egrep "Salmonella" /stor/scratch/Ochman/hassan/0318_AllTheBacteria/hq_dataset.species_calls.tsv | grep -v "Salmonella enterica" | cut -f1 | sed "s/^/fastANI -q \/stor\/work\/Ochman\/hassan\/protogene_extension\/Salmonella_list\/GCA_000210855.2_ASM21085v2_genomic.fna -r \/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\/\*\//g" | sed "s/$/\.fa/g" | awk -F '/' '{print $0" -o "$NF}' | sed "s/$/t/g" | sed "s/\.fat/_fastANI/g" | split -l 100
ls x* | sed "s/^/bash /g" > running.sh
./parallelize_run.sh

ls | egrep "fastANI$" | sed "s/^/cat /g" | sed "s/$/ >> Salmonella_all_fastANI/g" | bash #To prevent argument list too long error

awk -F '\t' '($2~"GCA_")' Salmonella_all_fastANI | cut -f2,3 | sed "s/\/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\///g" | sed 's/_/\t/2' | cut -f1,3 > Salmonella_accession_ANI.tsv
awk -F '\t' '($2!~"GCA_")' Salmonella_all_fastANI | cut -f2,3 | rev | cut -f1 -d '/' | rev | sed "s/\.fa\t/\t/g" >> Salmonella_accession_ANI.tsv
sort -k1 Salmonella_accession_ANI.tsv -o Salmonella_accession_ANI.tsv

#Add taxonomy:
grep "^GCA" Salmonella_accession_ANI.tsv | cut -f1 | grep -w -F -f - ../102824_GenBank_bacteria_assembly_summary.txt | sort -u | cut -f1,8 | sed "s/ /_/g" | sort -k1 | join -1 1 -2 1 - Salmonella_accession_ANI.tsv > Salmonella_taxonomy_accession_ANI.tsv
grep -v "^GCA" Salmonella_accession_ANI.tsv | cut -f1 | grep -w -F -f - /stor/scratch/Ochman/hassan/0318_AllTheBacteria/hq_dataset.species_calls.tsv | sed "s/ /_/g" | sort -k1 | join -1 1 -2 1 - Salmonella_accession_ANI.tsv >> Salmonella_taxonomy_accession_ANI.tsv
sort -nrk3 Salmonella_taxonomy_accession_ANI.tsv | sed "s/ /\t/g" > test && mv test Salmonella_taxonomy_accession_ANI.tsv

#Genus-excluded accessions:
awk -F '\t' '($3>89)' Salmonella_taxonomy_accession_ANI.tsv | grep "GCA_" | cut -f1 > Salmonella_genusexcluded_accesions.txt
awk -F '\t' '($3<89)' Salmonella_taxonomy_accession_ANI.tsv | egrep "Salmonella" | cut -f1 >> Salmonella_genusexcluded_accesions.txt

cd ../Salmonella_db
#Salmonella_excluded genomes:
ls ../genomes/ | grep "^GCA" | grep -F -v -f ../Salmonella_temp/Salmonella_genusexcluded_accesions.txt | sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" | sed "s/$/ >> Salmonella_excluded.genomes.faa/g" | bash
#Salmonella_excluded proteins:
ls ../proteins/ | grep "^GCA" | grep -F -v -f ../Salmonella_temp/Salmonella_genusexcluded_accesions.txt | sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/proteins\//g" | sed "s/$/ >> Salmonella_excluded.proteins.faa/g" | bash

#Salmonella_excluded_ORFs:
getorf -sequence Salmonella_excluded.genomes.faa -outseq Salmonella_excluded.genomes.getorf.all -table 1 -minsize 30 -find 3
seqkit fx2tab Salmonella_excluded.genomes.getorf.all | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/" > Salmonella_excluded.genomes.getorf.ATG_TTG_GTG
/stor/work/Ochman/hassan/tools/faTrans -stop Salmonella_excluded.genomes.getorf.ATG_TTG_GTG Salmonella_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa
pwd

#enterica_excluded genomes:
awk -F '\t' '($3<95.2&&$3>89)' ../Salmonella_temp/Salmonella_taxonomy_accession_ANI.tsv | grep -v "enterica.*enterica" | awk -F '\t' '($2!="Salmonella_enterica")' | grep "^GCA" | cut -f1 > Enterica_excluded_accessions.GBRS.txt
awk -F '\t' '($3<95.2&&$3>89)' ../Salmonella_temp/Salmonella_taxonomy_accession_ANI.tsv | grep -v "enterica.*enterica" | awk -F '\t' '($2!="Salmonella_enterica")' | grep -v "^GCA" | cut -f1 > Enterica_excluded_accessions.ATB.txt
sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" Enterica_excluded_accessions.GBRS.txt | sed "s/$/* >> Enterica_excluded.genomes.faa/g" | bash
sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\/\*\//g" Enterica_excluded_accessions.ATB.txt | sed "s/$/\.fa >> Enterica_excluded.genomes.faa/g" | bash

#enterica_excluded_proteins:
sed "s/^/ls \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" Enterica_excluded_accessions.GBRS.txt | sed "s/$/*/g" | bash > prodigal_interim
sed "s/^/ls \/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\/\*\//g" Enterica_excluded_accessions.ATB.txt | sed "s/$/\.fa/g" | bash >> prodigal_interim

#Run prodigal:
sed "s/^/prodigal -i /g" prodigal_interim | sed "s/$/ -f gff -o /g" | awk -F '/' '{print $0$NF}' | sed "s/ -f gff -o $//g" | sed "s/_genomic\.fna$/\.gff/g" | sed "s/\.fa$/\.gff/g" | split -l100
ls x* | sed "s/^/bash /g" > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh

#gff to gtf:
for i in $(ls *gff | rev | cut -f2- -d "." | rev); do grep -v "#" "$i".gff | awk -F '\t' '($3=="CDS")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,$1}' | rev | cut -f 2- -d "." | rev | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,"transcript_id \""$9"_"NR"\";gene_id \""$9"_"NR"\";"}' > "$i".gtf; done

#make space for the newly generated files:
mkdir gff
mkdir gtf
mkdir cds
mkdir proteins
mv *.gff gff
mv *.gtf gtf

#gtf to cds:
for i in $(ls gtf/GCA*gtf | rev | cut -f2- -d "." | rev | cut -f2- -d "/"); do cat gtf/"$i".gtf | gtf2bed | bedtools getfasta -s -name -fi /stor/scratch/Ochman/hassan/100724_Complete_Genomes/genomes/"$i"_genomic.fna -bed - > cds/"$i".cds.faa; done
for i in $(ls gtf/SA*gtf | rev | cut -f2- -d "." | rev | cut -f2- -d "/"); do cat gtf/"$i".gtf | gtf2bed | bedtools getfasta -s -name -fi /stor/scratch/Ochman/hassan/0318_AllTheBacteria/AllTheBacteria_OG/*/$i.fa -bed - > cds/$i.cds.faa; done

#cds to protein:
for i in $(ls cds/*faa | rev | cut -f3- -d '.' | rev | cut -f2- -d "/"); do /stor/work/Ochman/hassan/tools/faTrans -stop cds/$i.cds.faa proteins/$i.protein.faa; done
pwd

#Enterica_excluded_proteins:
cat proteins/*faa >> Enterica_excluded.proteins.faa

#Enterica_excluded_ORFs:
getorf -sequence Enterica_excluded.genomes.faa -outseq Enterica_excluded.genomes.getorf.all -table 1 -minsize 30 -find 3
seqkit fx2tab Enterica_excluded.genomes.getorf.all | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/" > Enterica_excluded.genomes.getorf.ATG_TTG_GTG
/stor/work/Ochman/hassan/tools/faTrans -stop Enterica_excluded.genomes.getorf.ATG_TTG_GTG Enterica_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa
pwd

makeblastdb -in Salmonella_excluded.genomes.faa -dbtype nucl -out Salmonella_excluded.genomes
makeblastdb -in Enterica_excluded.genomes.faa -dbtype nucl -out Enterica_excluded.genomes
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Salmonella_excluded.proteins.faa --db Salmonella_excluded.proteins
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Enterica_excluded.proteins.faa --db Enterica_excluded.proteins
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Salmonella_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa --db Salmonella_excluded.genomes.getorf.ATG_TTG_GTG.prot
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Enterica_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa --db Enterica_excluded.genomes.getorf.ATG_TTG_GTG.prot

#######################################################################################################################################
#Mycobacterium:
cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes/
mkdir Mycobacterium_temp && cd Mycobacterium_temp
cp /stor/work/Ochman/hassan/tools/parallelize_run.sh .
ls /stor/scratch/Ochman/hassan/100724_Complete_Genomes/genomes | grep "^GC" | cut -f1,2 -d "_" | sed "s/^/fastANI -q \/stor\/work\/Ochman\/hassan\/protogene_extension\/Mycobacterium_list\/H37Rv.fna -r \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" | sed "s/$/\*fna/g" | awk -F '/' '{print $0" -o "$NF}' | sed "s/$/t/g" | sed "s/\*fnat/_fastANI/g" | split -l 600
ls x* | sed "s/^/bash /g" > running.sh
./parallelize_run.sh
#For ATB: (could be made more non-redundant)
rm x*
egrep "Mycobacterium" /stor/scratch/Ochman/hassan/0318_AllTheBacteria/hq_dataset.species_calls.tsv | grep -v "Mycobacterium tuberculosis" | cut -f1 | sed "s/^/fastANI -q \/stor\/work\/Ochman\/hassan\/protogene_extension\/Mycobacterium_list\/H37Rv.fna -r \/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\/\*\//g" | sed "s/$/\.fa/g" | awk -F '/' '{print $0" -o "$NF}' | sed "s/$/t/g" | sed "s/\.fat/_fastANI/g" | split -l 100
ls x* | sed "s/^/bash /g" > running.sh
./parallelize_run.sh

ls | egrep "fastANI$" | sed "s/^/cat /g" | sed "s/$/ >> Mycobacterium_all_fastANI/g" | bash #To prevent argument list too long error

awk -F '\t' '($2~"GCA_")' Mycobacterium_all_fastANI | cut -f2,3 | sed "s/\/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\///g" | sed 's/_/\t/2' | cut -f1,3 > Mycobacterium_accession_ANI.tsv
awk -F '\t' '($2!~"GCA_")' Mycobacterium_all_fastANI | cut -f2,3 | rev | cut -f1 -d '/' | rev | sed "s/\.fa\t/\t/g" >> Mycobacterium_accession_ANI.tsv
sort -k1 Mycobacterium_accession_ANI.tsv -o Mycobacterium_accession_ANI.tsv

#Add taxonomy:
grep "^GCA" Mycobacterium_accession_ANI.tsv | cut -f1 | grep -w -F -f - ../102824_GenBank_bacteria_assembly_summary.txt | sort -u | cut -f1,8 | sed "s/ /_/g" | sort -k1 | join -1 1 -2 1 - Mycobacterium_accession_ANI.tsv > Mycobacterium_taxonomy_accession_ANI.tsv
grep -v "^GCA" Mycobacterium_accession_ANI.tsv | cut -f1 | grep -w -F -f - /stor/scratch/Ochman/hassan/0318_AllTheBacteria/hq_dataset.species_calls.tsv | sed "s/ /_/g" | sort -k1 | join -1 1 -2 1 - Mycobacterium_accession_ANI.tsv >> Mycobacterium_taxonomy_accession_ANI.tsv
sort -nrk3 Mycobacterium_taxonomy_accession_ANI.tsv | sed "s/ /\t/g" > test && mv test Mycobacterium_taxonomy_accession_ANI.tsv

#Mycobacterium-excluded accessions:

grep "Mycobacterium" Mycobacterium_taxonomy_accession_ANI.tsv | grep "GCA_" | cut -f1 > Mycobacterium_genusexcluded_accesions.txt
grep -v "Mycobacterium" Mycobacterium_taxonomy_accession_ANI.tsv | awk -F '\t' '($3>99)' | cut -f1 >> Mycobacterium_genusexcluded_accesions.txt

cd ../Mycobacterium_db
#Mycobacterium_excluded genomes:
ls ../genomes/ | grep "^GCA" | grep -F -v -f ../Mycobacterium_temp/Mycobacterium_genusexcluded_accesions.txt | sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" | sed "s/$/ >> Mycobacterium_excluded.genomes.faa/g" | bash

#Mycobacterium_excluded proteins:
ls ../proteins/ | grep "^GCA" | grep -F -v -f ../Mycobacterium_temp/Mycobacterium_genusexcluded_accesions.txt | sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/proteins\//g" | sed "s/$/ >> Mycobacterium_excluded.proteins.faa/g" | bash

#Mycobacterium_excluded_ORFs:
getorf -sequence Mycobacterium_excluded.genomes.faa -outseq Mycobacterium_excluded.genomes.getorf.all -table 1 -minsize 30 -find 3
seqkit fx2tab Mycobacterium_excluded.genomes.getorf.all | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/" > Mycobacterium_excluded.genomes.getorf.ATG_TTG_GTG
/stor/work/Ochman/hassan/tools/faTrans -stop Mycobacterium_excluded.genomes.getorf.ATG_TTG_GTG Mycobacterium_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa
pwd

#Tuberculosis_excluded genomes:
awk -F '\t' '($3<99&&$2~"Mycobacterium")' ../Mycobacterium_temp/Mycobacterium_taxonomy_accession_ANI.tsv | grep "^GCA" | cut -f1 > Tuberculosis_excluded_accessions.GBRS.txt
awk -F '\t' '($3<99&&$2~"Mycobacterium")' ../Mycobacterium_temp/Mycobacterium_taxonomy_accession_ANI.tsv | grep -v "^GCA" | cut -f1 > Tuberculosis_excluded_accessions.ATB.txt
sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" Tuberculosis_excluded_accessions.GBRS.txt | sed "s/$/* >> Tuberculosis_excluded.genomes.faa/g" | bash
sed "s/^/cat \/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\/\*\//g" Tuberculosis_excluded_accessions.ATB.txt | sed "s/$/\.fa >> Tuberculosis_excluded.genomes.faa/g" | bash

#Tuberculosis_excluded_proteins:
sed "s/^/ls \/stor\/scratch\/Ochman\/hassan\/100724_Complete_Genomes\/genomes\//g" Tuberculosis_excluded_accessions.GBRS.txt | sed "s/$/*/g" | bash > prodigal_interim
sed "s/^/ls \/stor\/scratch\/Ochman\/hassan\/0318_AllTheBacteria\/AllTheBacteria_OG\/\*\//g" Tuberculosis_excluded_accessions.ATB.txt | sed "s/$/\.fa/g" | bash >> prodigal_interim

#Run prodigal:
sed "s/^/prodigal -i /g" prodigal_interim | sed "s/$/ -f gff -o /g" | awk -F '/' '{print $0$NF}' | sed "s/ -f gff -o $//g" | sed "s/_genomic\.fna$/\.gff/g" | sed "s/\.fa$/\.gff/g" | split -l100
ls x* | sed "s/^/bash /g" > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh

#gff to gtf:
for i in $(ls *gff | rev | cut -f2- -d "." | rev); do grep -v "#" "$i".gff | awk -F '\t' '($3=="CDS")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,$1}' | rev | cut -f 2- -d "." | rev | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,"transcript_id \""$9"_"NR"\";gene_id \""$9"_"NR"\";"}' > "$i".gtf; done

#make space for the newly generated files:
mkdir gff
mkdir gtf
mkdir cds
mkdir proteins
mv *.gff gff
mv *.gtf gtf

#gtf to cds:
for i in $(ls gtf/GCA*gtf | rev | cut -f2- -d "." | rev | cut -f2- -d "/"); do cat gtf/"$i".gtf | gtf2bed | bedtools getfasta -s -name -fi /stor/scratch/Ochman/hassan/100724_Complete_Genomes/genomes/"$i"_genomic.fna -bed - > cds/"$i".cds.faa; done
for i in $(ls gtf/SA*gtf | rev | cut -f2- -d "." | rev | cut -f2- -d "/"); do cat gtf/"$i".gtf | gtf2bed | bedtools getfasta -s -name -fi /stor/scratch/Ochman/hassan/0318_AllTheBacteria/AllTheBacteria_OG/*/$i.fa -bed - > cds/$i.cds.faa; done

#cds to protein:
for i in $(ls cds/*faa | rev | cut -f3- -d '.' | rev | cut -f2- -d "/"); do /stor/work/Ochman/hassan/tools/faTrans -stop cds/$i.cds.faa proteins/$i.protein.faa; done
pwd

#Tuberculosis_excluded_proteins:
cat proteins/*faa >> Tuberculosis_excluded.proteins.faa

#Tuberculosis_excluded_ORFs:
getorf -sequence Tuberculosis_excluded.genomes.faa -outseq Tuberculosis_excluded.genomes.getorf.all -table 1 -minsize 30 -find 3
seqkit fx2tab Tuberculosis_excluded.genomes.getorf.all | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/" > Tuberculosis_excluded.genomes.getorf.ATG_TTG_GTG
/stor/work/Ochman/hassan/tools/faTrans -stop Tuberculosis_excluded.genomes.getorf.ATG_TTG_GTG Tuberculosis_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa
pwd

makeblastdb -in Mycobacterium_excluded.genomes.faa -dbtype nucl -out Mycobacterium_excluded.genomes
makeblastdb -in Tuberculosis_excluded.genomes.faa -dbtype nucl -out Tuberculosis_excluded.genomes
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Mycobacterium_excluded.proteins.faa --db Mycobacterium_excluded.proteins
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Tuberculosis_excluded.proteins.faa --db Tuberculosis_excluded.proteins
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Mycobacterium_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa --db Mycobacterium_excluded.genomes.getorf.ATG_TTG_GTG.prot
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Tuberculosis_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa --db Tuberculosis_excluded.genomes.getorf.ATG_TTG_GTG.prot
