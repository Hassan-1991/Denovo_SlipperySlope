# Dedup 99.99 → list of representative genomes; also write deduped 99.49 cluster file
awk '!seen[$1]++{print $2}' E.coli_pairwise_ANI.99.99.clusters > rep_genomes.txt
grep -F -f rep_genomes.txt E.coli_pairwise_ANI.99.49.clusters > E.coli_pairwise_ANI.99.49.clusters.dedup

# Clusters with >10 reps
cut -d' ' -f1 E.coli_pairwise_ANI.99.49.clusters.dedup | sort | uniq -c | awk '$1>9{print $2}' > clusters.txt

# Keep first 10 reps per retained cluster
for i in $(cat clusters.txt)
do
grep "^$i " E.coli_pairwise_ANI.99.49.clusters.dedup; done | awk '++seen[$1] <= 10' > nonred_genomes.withclusters.txt

#We now move these genomes to a new directory, place annotation files in a separate sub-directory, and give them informative names
#STARTING DIRECTORY:

/stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics

mkdir -p Ecoli_nonred_genomes
#mkdir -p Ecoli_nonred_gffs/genemarks2
#mkdir -p Ecoli_nonred_gffs/prodigal

awk '
{
  cluster=$1
  file=$2
  count[cluster]++
  newname="C"cluster"G"count[cluster]"_genomic.fna"
  src="../../Escherichia_coli_RefSeq_genomes/"file
  dst="Ecoli_nonred_genomes/"newname
  print "cp \""src"\" \""dst"\""
}
' nonred_genomes.withclusters.txt | bash

for f in *.fna; do
  prefix=$(basename "$f" | cut -f1 -d "_")
  awk -v pfx="$prefix" '
    BEGIN{n=0}
    /^>/{
      if (tolower($0) ~ /plasmid/) { n++; print ">" pfx "_plasmid_" n }
      else                        { print ">" pfx "_chromosome" }
      next
    }
    {print}
  ' "$f" > "${f%.fna}.renamed.fna"
done

for i in $(ls *_genomic.renamed.fna | cut -f1 -d "_")
do
echo "prodigal -i ${i}_genomic.renamed.fna -f gff -o ${i}_prodigal.gff"
done | split -l10

ls x* | sed "s/^/bash /g" > prodigal_run
/stor/work/Ochman/hassan/tools/parallelize_run.sh prodigal_run

#Genemarks2:

for i in $(ls *_genomic.renamed.fna | cut -f1 -d "_")
do
dir="${i}_gms2"
echo "mkdir -p \"$dir\" && cd \"$dir\" && /stor/work/Ochman/hassan/tools/gms2_linux_64/gms2.pl --seq \"../${i}_genomic.renamed.fna\" --genome-type bacteria --gcode 11 --output \"${i}_genemarks2.gff\" --format gff3"
done > genemarks_run

/stor/work/Ochman/hassan/tools/parallelize_run.sh genemarks_run

#Move to appropriate directories
