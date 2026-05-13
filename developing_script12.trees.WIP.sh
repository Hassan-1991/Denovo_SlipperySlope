#These are steps that need to be incorporated way earlier in the pipeline

#Make a tree with the clusters
mkdir -p phame
awk '!seen[$1]++' nonred_genomes.withclusters.txt | sort -k2 > temp
sort -k1 renames_tally.tsv | join -1 1 -2 2 - temp | cut -f2 -d " " | sed "s/^/cp Ecoli_nonred_genomes\//g" | sed "s/$/ phame/g" | sed "s/genomic/genomic.renamed/g" | bash

#ALL:
sort -k2 nonred_genomes.withclusters.txt | join -1 2 -2 1 - renames_tally.tsv | cut -f2- -d " " | awk '{print "cp Ecoli_nonred_genomes\/"$2,"phame_all/"$1"_"$2}' | sed "s/genomic/genomic.renamed/g" | bash

#Download an outgroup (E. fergusonii):
cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes
grep "Escherichia fergusonii" 100225_RefSeq_assembly_summary.txt | grep "Complete Genome" | grep -i "reference genome" | cut -f20 | awk -F '/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed "s/^/wget /g" | bash
gunzip *_genomic.fna.gz

cp GCF_020097475.1_ASM2009747v1_genomic.fna Ecoli_genomics/Ecoli.clusters.0.1.gaps/phame/

cd Ecoli_genomics/Ecoli.clusters.0.1.gaps/phame/

mkdir -p refdir
mkdir -p workdir

mv *fna refdir

conda activate phame
phame Ecoli.ctl

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

#Rename trees:
awk -F'\t' 'NR==FNR{m[$1]=$3; next} {for(k in m) gsub(k"_", m[k]"_"k"_"); gsub(/_genomic_renamed/, ""); print}' genome_cluster.contigs.tsv phame/workdir/results/trees/RAxML_bestTree.Ecoli_all



#Just for E. coli, since we have the poppunk info already, let's attach that:

#First, run pairwise fastANI between the two sets of genomes:
#!/bin/bash

qdir="/stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/Ecoli_nonred_genomes"
rdir="/stor/work/Ochman/hassan/MS_Ecoli_ORFans_Ch3/Ecoli_genomes"
outdir="/stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/two_genome_set_comparisons"

mkdir -p "$outdir"

: > fastANI_two_sets_mastercode.sh

for q in "$qdir"/*fna; do
    qbase=$(basename "$q")
    qbase=${qbase%.*}

    for r in "$rdir"/*fasta; do
        rbase=$(basename "$r")
        rbase=${rbase%.*}

        echo "fastANI -q \"$q\" -r \"$r\" -o \"$outdir/${qbase}_${rbase}.tsv\"" >> fastANI_two_sets_mastercode.sh
    done
done

#Huge time investment. Once that's run in a parallelized way, let's pick the sequence in the poppunk set that the RefSeq set has the best hit with:

qdir="/stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/Ecoli_nonred_genomes"
tdir="/stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/two_genome_set_comparisons"

mkdir -p best_per_query

for q in "$qdir"/*_genomic.renamed.fna; do
    qbase=$(basename "$q" .fna)

    find "$tdir" -name "${qbase}_*.tsv" -print0 |
    xargs -0 awk -F '\t' '
        $3 > best {best=$3; line=$0}
        END {if (line) print line}
    ' > "best_per_query/${qbase}.best.tsv" &
done

wait

#Extract phylogroup and lineage info for the popppunk genomes, attach:

cd best_per_query

cat *best.tsv |
sed "s/\/stor\/work\/Ochman\/hassan\/MS_Ecoli_ORFans_Ch3\/Ecoli_genomes\///g" |
sed "s/.fasta//g" | sort -k2 |
join -1 2 -2 1 - 500_ipp_lineagedesignations.sorted.tsv |
awk '{print $2,$6"@"$NF}' | rev | cut -f1 -d "/" | rev | sed "s/_genomic.renamed.fna//g" | sort -k1 |
join -1 1 -2 1 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/genome_cluster.contigs.tsv | sed "s/ /\t/g" > /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/genome_phylogroup_contig_lineage.tsv

#Rename the tree:

sed "s/_genomic_renamed//g" phame/workdir/results/trees/RAxML_bestTree.Ecoli_all | perl -F'\t' -lane '
BEGIN {
    open M, "genome_phylogroup_contig_lineage.tsv" or die $!;
    while (<M>) {
        chomp;
        @x = split /\t/;
        $map{$x[0]} = "$x[0]_$x[1]_$x[3]";
    }
}
for $k (keys %map) {
    $q = quotemeta($k);
    s/([\(,])$q(?=:)/$1$map{$k}/g;
}
print;
' > Ecoli_RaxML_tree.nwk

#Pretty much it. Rotate branches as needed, jot down lineage order. Khalas
