#Download all E. coli complete genomes from RefSeq:

# 1) get RefSeq assembly summary
wget -q -O 100225_RefSeq_assembly_summary.txt https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

# 2) build URL list (Complete Genome + organism contains "Escherichia coli")
awk -F '\t' '($8 ~ /Escherichia coli/ && $12=="Complete Genome"){print $20}' 100225_RefSeq_assembly_summary.txt |
awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' > 100225_RefSeq_Ecoli_genome_URLs.txt

# 3) download into a dedicated folder
mkdir -p Escherichia_coli_RefSeq_genomes
mv 100225_RefSeq_Ecoli_genome_URLs.txt Escherichia_coli_RefSeq_genomes/
cd Escherichia_coli_RefSeq_genomes

wget -q -i 100225_RefSeq_Ecoli_genome_URLs.txt -c
gunzip -f *.gz

cd ..

#fastANI

# 1) Make a mastercode of all-vs-all fastANI comparisons

#List genomes 1-each-line, put it in an array called genomes
mapfile -t genomes < <(ls -1 *genomic.fna)

: > fastANI_allvall_mastercode.sh #Create file

n=0
for q in "${genomes[@]}"; do
  for r in "${genomes[@]}"; do
    [[ "$q" == "$r" ]] && continue #skip if files are same
    n=$((n+1))
    out="fastani_${n}.tsv" #Each line's output directs to a uniquely numbered output file
    echo "fastANI -q \"$q\" -r \"$r\" -o \"$out\"" >> fastANI_allvall_mastercode.sh
  done
done

# 2) Parallelize this across 176 threads in two servers

total_chunks=176
lines=$(wc -l < fastANI_allvall_mastercode.sh)
per=$(( (lines + total_chunks - 1) / total_chunks ))   #Calculate number of lines that go into per file

split -d -a 3 -l "$per" fastANI_allvall_mastercode.sh fastani_chunk_ #Split the mastecode file into 176 smaller code files

ls -1 fastani_chunk_* | head -104 | awk '{print "bash "$0}' > ochmcomp02_running.sh 
ls -1 fastani_chunk_* | tail -72  | awk '{print "bash "$0}' > ochmcomp01_running.sh #Two files to run, given 2 servers

#Helper script for parallelization:

#!/usr/bin/env bash
set -euo pipefail

runfile="${1:?usage: parallelize_run.sh RUNFILE}"

while IFS= read -r line; do
  screen -dmS "job_$RANDOM" bash -lc "$line"
done < "$runfile"

./parallelize_run.sh ochmcomp02_running.sh
./parallelize_run.sh ochmcomp01_running.sh #Actual run

#Compile all of these results into one file:

find . -maxdepth 1 -name 'fastani_chunk_*' | xargs cat > all_fastANI_comparisons.tsv
