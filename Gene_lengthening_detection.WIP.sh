#Starting from developing_script6.blastp_parsing:
#all_proteins_vs_GBRS_annotated.tsv
#And /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_extragenus_protein_contig_taxa.reduced.noescherichia.tsv from preparing_databases in the other repository:

#Assign taxa to all hits in the Ecoli vs all GBRS annotated diamond blast file:
sort -k2 all_proteins_vs_GBRS_annotated.tsv | join -1 2 -2 2 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_extragenus_protein_contig_taxa.reduced.noescherichia.tsv > all_proteins_vs_GBRS_annotated.taxa.tsv

#Identify small proteins (<101AA) which have at least 1 hit against a protein >250AA in length:
sed "s/ /\t/g" all_proteins_vs_GBRS_annotated.taxa.tsv |
awk -F '\t' '($14<101)' | #Small protein
awk '($16<0.001)' | awk '($15>250)' | #at least one hit against a protein >250AA
awk -F '\t' '($5>80)' | #Virtually the entire length of the small protein must be a part of the longer protein
cut -f2 | sort -u > Ecoli_RefSeq.proteins_with_longerhits.txt

#In parallel, identify long proteins (>250AA) which have at least 1 hit against a protein <100AA in length:
seqkit fx2tab /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/Ecoli_nonred_gffs/reducing_redundancies_expt/Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa |
awk '{print $1,length($2)}' | awk '($2>250)' | #Select long proteins
cut -f1 -d " " | sort -u |
grep -w -f - /stor/work/Ochman/hassan/MS_denovoreview_Ch5/all_proteins_vs_GBRS_annotated.taxa.tsv | #Identify corresponding blast lines
awk '(($6/$15)>0.8)' | #Alignment length/subject length > 0.8——same logic, entire small protein must be found within large protein
awk '($15<100)' > all_proteins_vs_GBRS_annotated.genelengtheningcandidates.tsv #at least one hit <100AA

cut -f2 -d " " all_proteins_vs_GBRS_annotated.genelengtheningcandidates.tsv | sort | uniq -c | awk '($1!=1)' | #Get rid of cases in which the short length appears just once——annotation error or de novo genome-specific mutation
sort -u | awk '{print $NF}' > Ecoli_RefSeq.proteins_with_shorterhits.txt

#Inside new directory
#Compile all homologous sequences into one file and reduce sequence redundancy with usearch:

for i in $(cat ../Ecoli_RefSeq.proteins_with_longerhits.txt)
do
echo 'grep -w "'"${i}"'" ../all_proteins_vs_GBRS_annotated.taxa.tsv | sed "s/ /\t/g" | cut -f1 | sort -u > '"${i}"'.hits.txt' #Get all homolog identifiers
echo 'faSomeRecords /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.proteins.faa '"${i}"'.hits.txt '"${i}"'.hits.faa' #Get their sequences. Preparing databses from other repo
echo 'usearch -sortbylength '"${i}"'.hits.faa -fastaout '"${i}"'.hits.sorted.faa -minseqlength 1' #sort in order to cluster
echo 'usearch -cluster_smallmem '"${i}"'.hits.sorted.faa -id 0.9 -centroids '"${i}"'.hits.centroids.faa -uc '"${i}"'.hits.uc' #Cluster at 90% pid, entire length of protein
echo 'seqkit fx2tab '"${i}"'.hits.centroids.faa | sed "s/\t$//g" | awk '\''{print ">"$1"\t"$NF}'\'' | sed "s/\t/\n/g" > '"${i}"'.mafft.input.faa' #Linearize
#echo 'grep -w -A1 "'"${i}"'" /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/Ecoli_nonred_gffs/reducing_redundancies_expt/Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa >> '"${i}"'.mafft.input.faa' #Add in query/gocal sequence
done

#For longer hits, same logic:

for i in $(cat ../Ecoli_RefSeq.proteins_with_shorterhits.txt)
do
echo 'grep -w "'"${i}"'" ../all_proteins_vs_GBRS_annotated.taxa.tsv | sed "s/ /\t/g" | cut -f1 | sort -u > '"${i}"'.hits.txt'
echo 'faSomeRecords /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.proteins.faa '"${i}"'.hits.txt '"${i}"'.hits.faa'
echo 'usearch -sortbylength '"${i}"'.hits.faa -fastaout '"${i}"'.hits.sorted.faa -minseqlength 1'
echo 'usearch -cluster_smallmem '"${i}"'.hits.sorted.faa -id 0.9 -centroids '"${i}"'.hits.centroids.faa -uc '"${i}"'.hits.uc'
echo 'seqkit fx2tab '"${i}"'.hits.centroids.faa | sed "s/\t$//g" | awk '\''{print ">"$1"\t"$NF}'\'' | sed "s/\t/\n/g" > '"${i}"'.mafft.input.faa'
done

#Now the controversial step: In a phylogeny-blind way, count up number of short vs. long homologous sequences
#If the number of short homologs exceed that of long homologs, tentative fusion candidate, retain.
#In order to not make this filter very strict, we should try to minimize the number of long homologs
#As long as #of <100AA homologs exceed the number of #250 homologs, they can be considered ancestral

for f in $(sed "s/$/.mafft.input.faa/g" Ecoli_RefSeq.proteins_with_shorterhits.250AA.txt)
do
    i=${f%.mafft.input.faa}

    seqkit fx2tab "$f" |
    awk -v id="$i" '
    {
        len = length($2)
        if (len < 100) short++
        else if (len > 250) long++
    }
    END {
        if (short > long) print id
    }'
done > extragenus_candidates.short_in_Ecoli.txt

#Likewise for the proteins with long E. coli homologs:
for f in $(sed "s/$/.mafft.input.faa/g" Ecoli_RefSeq.proteins_with_shorterhits.txt)
do
    i=${f%.mafft.input.faa}

    seqkit fx2tab "$f" |
    awk -v id="$i" '
    {
        len = length($2)
        if (len < 100) short++
        else if (len > 250) long++
    }
    END {
        if (short > long) print id
    }'
done > extragenus_candidates.long_in_Ecoli.txt

#Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa comes from developing_script11.reducingredundancy
#Make an interim file with homolog ID, sequence, focal sequence protein, homolog genome, genus

for i in $(cat extragenus_candidates.short_in_Ecoli.txt)
do
echo 'awk '\''NR==FNR{qlen=length($2); next} {print $1,$NF,qlen}'\'' <(grep -A1 "'"${i}"'" /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/Ecoli_nonred_gffs/reducing_redundancies_expt/Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa | seqkit fx2tab) <(seqkit fx2tab '"${i}"'.hits.faa) | sort -k1 | join -1 1 -2 2 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_extragenus_protein_contig_taxa.reduced.noescherichia.tsv > '"${i}"'.hits.length.taxa.tsv'
done | split -l1 #to parallelize this time-consuming step

#Now for the complicated step
#First, transform the file we just made to make a three-column file which lists how many times proteins of each length appear in each genus
#Genus - length - number of times that lengthed protein appears in that genus
#This gives us the *length_frequency.tsv files

#In parallel, we also make a list of all the genera that the homologs for each gene belong to
#Using a python script, we extract the tree for these genera from the gtdb newick file
#The tree is rooted according to the genus that's furthest removed from Escherichia in the collection we're working with
#This gives us the *gtdb.rooted_by_farthest_from_Escherichia.tree files

#Using a R script now, we label the leaves of the tree according to how many short/long genes are present in the genus
#If all genomes in the genus are short or long, they get colored blue or red, respectively
#If there's a mix, they're colored purple, with a label stating what fraction of genomes in that genus are in the minority
#Escherichia (being the focal genus) gets labeled black

for f in *hits.length.taxa.tsv
do
  prefix=${f%.hits.length.taxa.tsv}

  cat > "${prefix}.plot_pipeline.sh" <<EOF

#!/usr/bin/env bash
set -euo pipefail

awk '{print \$NF}' ${prefix}.hits.length.taxa.tsv | sort -u | sed "1s/^/Escherichia\\n/g" > ${prefix}.genuslist.txt

focal_len=\$(awk 'NR==1{print \$3; exit}' ${prefix}.hits.length.taxa.tsv)

cut -f2,5 -d " " ${prefix}.hits.length.taxa.tsv |
awk '{print \$NF,length(\$1)}' |
sort | uniq -c |
awk '{print \$2,\$3,\$1}' |
sort -nrk2 |
sed "s/ /\\t/g" |
sed "1s/^/Escherichia\\tfocal\\t\${focal_len}\\n/g" \
> ${prefix}.length_frequency.tsv

python make_gtdb_genus_tree.py ${prefix}.genuslist.txt

Rscript plot_length_tree.R ${prefix}
EOF

  chmod +x "${prefix}.plot_pipeline.sh"
done


