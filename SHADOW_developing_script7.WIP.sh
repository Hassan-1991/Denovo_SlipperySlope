#Regular blastn and parse

#Regular blastn, without regard to flanks:

#Re: options:
#The num_alignments and num_descriptions reflect the number of target sequences (contigs) in each database
#We set -outfmt to 0, as we'll be processing the results with mview later
#evalue is kept at an arbitrarily high number——this search is not restricted by evalue cutoff

#We make the pangenome blast db:
makeblastdb -in Ecoli.genome.pangenomedb.fna -dbtype nucl -out Ecoli.genome.pangenomedb

mkdir shadow_analysis
cp Ecoli.genome.pangenomedb* shadow_analysis
cp Ecoli_queries.CDS.fna shadow_analysis
cd shadow_analysis

#Escherichia-excluded:
blastn -query Ecoli_queries.CDS.fna -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.genomes -outfmt 0 -num_threads 72 -num_descriptions 109734 -num_alignments 109734 -evalue 200000 -out Ecoli_extragenus_regular_blastn
#Ecoli-excluded:
blastn -query Ecoli_queries.CDS.fna -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Ecoli_excluded.genomes -outfmt 0 -num_threads 72 -num_descriptions 195119 -num_alignments 195119 -evalue 200000 -out Ecoli_intragenus_regular_blastn
#Pangenome:
blastn -query Ecoli_queries.CDS.fna -db Ecoli.genome.pangenomedb -outfmt 0 -num_threads 72 -num_descriptions 1761 -num_alignments 1761 -evalue 200000 -out Ecoli_pangenome_regular_blastn

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

find . -name "*_blastn" -print0 | xargs -0 grep "No hits found" | grep -v "regular" | cut -f1 -d ":" | sed "s/^/rm /g" | bash

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

