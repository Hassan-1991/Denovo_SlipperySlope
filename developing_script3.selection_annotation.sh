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

#Save the mappings between old and new names for later:

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
' | sed "s/ /\//g" | cut -f5,7 -d "/" | sed "s/\//\t/g" | tr -d "\"" > renames_tally.tsv

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

#Make a tree with the clusters
mkdir -p phame
awk '!seen[$1]++' nonred_genomes.withclusters.txt | sort -k2 > temp
sort -k1 renames_tally.tsv | join -1 1 -2 2 - temp | cut -f2 -d " " | sed "s/^/cp Ecoli_nonred_genomes\//g" | sed "s/$/ phame/g" | sed "s/genomic/genomic.renamed/g" | bash

#Download an outgroup (E. fergusonii):
cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes
grep "Escherichia fergusonii" 100225_RefSeq_assembly_summary.txt | grep "Complete Genome" | grep -i "reference genome" | cut -f20 | awk -F '/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed "s/^/wget /g" | bash
gunzip *_genomic.fna.gz

cp GCF_020097475.1_ASM2009747v1_genomic.fna Ecoli_genomics/Ecoli.clusters.0.1.gaps/phame/

cd Ecoli_genomics/Ecoli.clusters.0.1.gaps/phame/

mkdir -p refdir
mkdir -p workdir

mv *fna refdir

#Edit .ctl file appropriately:

   refdir = /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/phame/refdir/  # directory where reference (Complete) files are located
  workdir = /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/phame/workdir/ # directory where contigs/reads files are located and output is stored

reference = 2  # 0:pick a random reference from refdir; 1:use given reference; 2: use ANI based reference
  reffile =  # reference filename when option 1 is chosen

  project = Ecoli  # main alignment file name

  cdsSNPS = 0  # 0:no cds SNPS; 1:cds SNPs, divides SNPs into coding and non-coding sequences, gff file is required

  buildSNPdb = 0 # 0: only align to reference 1: build SNP database of all complete genomes from refdir

FirstTime = 1  # 1:yes; 2:update existing SNP alignment, only works when buildSNPdb is used first time to build DB

     data = 0  # *See below 0:only complete(F); 1:only contig(C); 2:only reads(R);
               # 3:combination F+C; 4:combination F+R; 5:combination C+R;
               # 6:combination F+C+R; 7:realignment  *See below
    reads =  # 1: single reads; 2: paired reads; 3: both types present;

     tree = 1  # 0:no tree; 1:use FastTree; 2:use RAxML; 3:use both;
bootstrap = 1  # 0:no; 1:yes;  # Run bootstrapping  *See below
        N = 100  # Number of bootstraps to run *See below

PosSelect = 1  # 0:No; 1:use PAML; 2:use HyPhy; 3:use both # these analysis need gff file to parse genomes to genes

     code = 0  # 0:Bacteria; 1:Virus; 2: Eukarya # Bacteria and Virus sets ploidy to haploid

    clean = 1  # 0:no clean; 1:clean # remove intermediate and temp files after analysis is complete

  threads = 72  # Number of threads to use

   cutoff = 0.1  # Linear alignment (LA) coverage against reference - ignores SNPs from organism that have lower cutoff.

   * When using data option 1,2,5 need a complete reference to align/map to.
   * Use data option 7 when need to extract SNPs using a sublist of already aligned genomes.

