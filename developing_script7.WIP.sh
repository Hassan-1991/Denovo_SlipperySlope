#Regular blastn and parse

#Regular blastn, without regard to flanks:

#Re: options:
#The num_alignments and num_descriptions reflect the number of target sequences (contigs) in each database
#We set -outfmt to 0, as we'll be processing the results with mview later
#evalue is kept at an arbitrarily high number——this search is not restricted by evalue cutoff

#We make the pangenome blast db:
makeblastdb -in Ecoli.genome.pangenomedb.fna -dbtype nucl -out Ecoli.genome.pangenomedb

#Escherichia-excluded:
blastn -query step1_ORFans.CDS.fna -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.genomes -outfmt 0 -num_threads 72 -num_descriptions 109734 -num_alignments 109734 -evalue 200000 -out Ecoli_extragenus_regular_blastn
#Ecoli-excluded:
blastn -query step1_ORFans.CDS.fna -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Ecoli_excluded.genomes -outfmt 0 -num_threads 72 -num_descriptions 195119 -num_alignments 195119 -evalue 200000 -out Ecoli_intragenus_regular_blastn
#Pangenome:
blastn -query step1_ORFans.CDS.fna -db Ecoli.genome.pangenomedb -outfmt 0 -num_threads 72 -num_descriptions 1761 -num_alignments 1761 -evalue 200000 -out Ecoli_pangenome_regular_blastn

mkdir regular
mv *_regular_blastn regular
cd regular

#We split each of the blast files to generate query-specific files
#Each resulting file is named after the query sequence

for i in extragenus intragenus pangenome
do
awk -v prefix="${i}_" '
/^Query=/ {
    close(file)
    q = $0
    sub(/^Query= /, "", q)
    file = prefix q "_blastn"
    next
}
{
    if (file) print > file
}' Ecoli_"$i"_regular_blastn
done

#Delete the query-specific files that contain "No hits found", i.e. no blast alignment

grep "No hits found" *_blastn | grep -v "regular" | cut -f1 -d ":" | sed "s/^/rm /g" | bash

#For subsequent processing with mview, we then prepend and append blast file fluff to beginning and end of each file

for i in extragenus intragenus pangenome
do
header=$(head -14 Ecoli_"$i"_regular_blastn) #assign header to variable
footer=$(tail -11 Ecoli_"$i"_regular_blastn) #assign footer to variable
for file in $(ls "$i"_*blastn); do
    # Prepend the header to the file
    { echo "$header"; cat "$file"; } > temp_file && mv temp_file "$file"

    # Append the footer to the file
    { cat "$file"; echo "$footer"; } > temp_file && mv temp_file "$file"
done
done

#LET'S PAUSE HERE

#Now we convert the blast output to a parseable mview file
export PERL5LIB=/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/lib/ #necessary to run this first for mviewed to run

#Since mview takes a while, we put the code in "echo" to parallelize it
for i in $(ls *_blastn | grep -v "regular" | rev | cut -f2- -d "_" | rev)
do
echo "/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/bin/mview -in blast ${i}_blastn > ${i}_blastn_mviewed"
done > running.sh

/stor/work/Ochman/hassan/tools/parallelize_run.sh running.sh #helper script used for parallelization

#We now parse the mviewed file according to some criteria
#We get rid of any case where the total combined length of all alignments doesn't cover 50% of query length
#And also targets whose alignments don't cover the query by at least 50%
#Since there's no e-value cutoff, it's necessary to utilize a query coverage cutoff to get rid of spurious hits
#This might miss out on very small regions of homologous non-coding sequences
#But that's fine, some of this information is recovered at the next stage of flank/synteny-based analysis

for i in $(ls *mviewed | grep -v "regular" | rev | cut -f3- -d "_" | rev)
do
querylength=$(echo $i | cut -f2- -d "_" | grep -w -F -A1 -f - ../step1_ORFans.CDS.fna | grep -v "^>" | awk '{print length($0)}') #length of query
ratio=$(tail -n+8 "$i"_blastn_mviewed | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}') #ratio of the total combined length of any and all alignments and the query length
if (( $(echo "$ratio > 0.5" | bc -l) )) #Only consider the total length of the query covered by any and all alignments is over 50% of the query length
then
tail -n+9 "$i"_blastn_mviewed | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2>50)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" > "$i"_blastn_seq.faa #Massaging
fi
done

#We then put the different massaged blast sequences in the same file:
for i in $(ls *_blastn_seq.faa | cut -f2- -d "_" | rev | cut -f3- -d "_" | rev)
do
cat extragenus_"$i"_blastn_seq.faa intragenus_"$i"_blastn_seq.faa pangenome_"$i"_blastn_seq.faa >> "$i"_blastn.seq.faa
done

#DO we even need the next part? Revisit.

####################
####################
####################
####################
####################


#In the next step, we remove all cases where there was a hit against a protein

#Attach taxonomy info to each file
#Make an interim file that contains the taxonomy information for all target contigs
#Otherwise searching the massive database would take a long time
cat *_blastn.seq.faa | grep "^>" | tr -d ">" | sort -u > alltaxa.compiled_intervalinfo.txt #intervalinfo is a vestigial name but we roll with it
cat /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_extragenus_genome_contig_taxa.reduced.noescherichia.tsv /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_intragenus_genome_contig_taxa.tsv /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/genome_cluster.contigs.tsv | grep -w -F -f alltaxa.compiled_intervalinfo.txt - | cut -f2- > alltaxa.compiled_intervalinfo.interim
sort -k1 alltaxa.compiled_intervalinfo.interim -o alltaxa.compiled_intervalinfo.interim

#Attach taxonomy information from the nascent file to a linearized blast alignment file
#The script is wrapped in echo to parallelize
for i in $(ls *_blastn.seq.faa | rev | cut -f2- -d "_" | rev)
do
echo "seqkit fx2tab ${i}_blastn.seq.faa | sed \"s/\t$//g\" | sort -k1 > ${i}_blastn.seq.tsv" #linearize blast sequence file
echo "cut -f1 ${i}_blastn.seq.tsv | sort -u > ${i}.temp" #get target contig names
echo "grep -w -F -f ${i}.temp alltaxa.compiled_intervalinfo.interim | sort -k1 | join -1 1 -2 1 - ${i}_blastn.seq.tsv > ${i}_blastn.seq.interim.tsv" #add taxonomy info to each target sequence
done | split -l3

ls x* | sed "s/^/bash /g" > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh running.sh #helper script used for parallelization

#If a target had a protein-based hit, remove it from the list of targets
for i in $(ls *_blastn.seq.interim.tsv | rev | cut -f2- -d "_" | rev)
do
echo "$i" | sed "s/$/(/g" | grep -F -f - ../all_genes_of_interest.presence.tsv | cut -f2 -d " " | grep -v -w -F -f - "$i"_blastn.seq.interim.tsv |
awk '{print ">regular_"$2"_"$1,$3}' | sed "s/ /\n/g" > "$i"_blastn.seq.regular.faa #mark the sequences "regular" to indicate how these hits were found
done

wc -l *_blastn.seq.regular.faa | awk '($1==0)' | grep -v "total" | rev | cut -f1 -d " " | rev | sed "s/^/rm /g" | bash #If any homology files are empty after the removal step, get rid of them
