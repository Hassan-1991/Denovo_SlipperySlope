mkdir flank
cp Ecoli.nr.filtered.prot.clusters.longestrepresentative.silixcollapsed.tsv flank
cd flank

#The next bit of script progresses in a few steps
#1. For each member in a gene family of interest, get its upstream and downstream gene
#2. Make a five-column tab-separated that lists the gene of interest, upstream gene, its cluster, downstream gene, its cluster
#3. Concatenate these rows, each corresponding to a member of the gene familiy, into a single file
#4. If any row contains fewer than five columns, that means that gene does not have both upstream and downstream flanking genes
#5. Retain each unique pair of family name for the upstream and downstream genes

#For each gene in the list of genes of interest (e.g., ORFan candidates)
#We loop over each of the members in their family
#Since the code is time-consuming, I put it in an "echo" wrapper below
#I'm showing an example for one gene of interest to demonstrate how the code looks
for j in $(cat ../step1_ORFans.txt)
do
    echo "for i in \$(awk '(\$1==\"$j\")' Ecoli.nr.filtered.prot.clusters.longestrepresentative.silixcollapsed.tsv | cut -f2)" >> "${j}.sh"
    echo "do" >> "${j}.sh"
    echo "    genome_name=\$(echo \$i | cut -f2 -d \"@\" | cut -f1 -d \"_\")" >> "${j}.sh"
    echo "    gene_ortholog_gtf=\$(grep -w \$i ../nr_gtf/\$genome_name.filtered.nr.gtf)" >> "${j}.sh"
    echo "    echo \$i > ${j}_temp" >> "${j}.sh"
    echo "    echo \$gene_ortholog_gtf | sed \"s/ /%/9\" | sed \"s/ /%/9\" | sed \"s/ /\\t/g\" | sed \"s/%/ /g\" | bedtools closest -a - -b ../nr_gtf/\$genome_name.filtered.nr.gtf -D a -id -io | awk -F '\\t' '(\$12==\"CDS\")' | cut -f 6 -d \"\\\"\" >> \"\$i\"_upstream" >> "${j}.sh"
    echo "    cat \"\$i\"_upstream >> ${j}_temp" >> "${j}.sh"
    echo "    grep -w -F -f \"\$i\"_upstream Ecoli.nr.filtered.prot.clusters.longestrepresentative.silixcollapsed.tsv | cut -f1 | sort -u >> ${j}_temp" >> "${j}.sh"
    echo "    rm \"\$i\"_upstream" >> "${j}.sh"
    echo "    echo \$gene_ortholog_gtf | sed \"s/ /%/9\" | sed \"s/ /%/9\" | sed \"s/ /\\t/g\" | sed \"s/%/ /g\" | bedtools closest -a - -b ../nr_gtf/\$genome_name.filtered.nr.gtf -D a -iu -io | awk -F '\\t' '(\$12==\"CDS\")' | cut -f 6 -d \"\\\"\" >> \"\$i\"_downstream" >> "${j}.sh"
    echo "    cat \"\$i\"_downstream >> ${j}_temp" >> "${j}.sh"
    echo "    grep -w -F -f \"\$i\"_downstream Ecoli.nr.filtered.prot.clusters.longestrepresentative.silixcollapsed.tsv | cut -f1 | sort -u >> ${j}_temp" >> "${j}.sh"
    echo "    rm \"\$i\"_downstream" >> "${j}.sh"
    echo "    if [ \"\$(wc -l < ${j}_temp)\" -eq 5 ]; then" >> "${j}.sh"
    echo "        paste -sd, ${j}_temp | sed \"s/,/\\t/g\" >> ${j}_flank_analysis" >> "${j}.sh"
    echo "        awk -F '\\t' 'NF==5' ${j}_flank_analysis | cut -f3,5 | sort -u > ${j}_flanks" >> "${j}.sh"
    echo "    fi" >> "${j}.sh"
    echo "done" >> "${j}.sh"
done

ls *.sh | sed "s/^/bash /g" > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh running.sh

#Follow next steps. Maybe don't search for homologs, instead put the regular and flanked regions together and search for homologs
#By, y'know, extracting ORFs and searching

#Make a helper file concatenating the sequences of all flanking genes so next steps go quicker
cat *_flanks | sed "s/\t/\n/g" | sort -u | grep -A1 --no-group-separator -F -f - ../Ecoli.nr.filtered.CDS.faa > all_flanks.CDS.faa

#Geneflanks:

#Extracting flanking genes is straightforward
#For each unique flanking gene pair, assign a line number
#Extract their CDSs
#The identifier should be of the form: GeneOfInterest_linenumber_up/down_FlankingGeneFamily

for j in $(ls *_flanks | rev | cut -f2- -d "_" | rev)
do
for i in $(cat -n "$j"_flanks | sed "s/^ *//g" | sed "s/\t/,/g")
do
linenumber=$(echo $i | cut -f1 -d ',')
echo $i | cut -f2 -d ',' | grep --no-group-separator -A1 -f - all_flanks.CDS.faa | sed "s/>/>"$j"_"$linenumber"_up_/g" >> geneflanks.faa
echo $i | cut -f3 -d ',' | grep --no-group-separator -A1 -f - all_flanks.CDS.faa | sed "s/>/>"$j"_"$linenumber"_down_/g" >> geneflanks.faa
done
done

#Proxflanks:

#To get sequences immediately upstream and downstream of the gene of interest:
#Identify a genome——any genome——just pick one at random——where the gene has both flanking genes
#Extract +/-500bp

for j in $(ls *_flanks | rev | cut -f2- -d "_" | rev)
do
#The next line takes one ortholog at random from each unique pair of flanking genes
for i in $(awk -F '\t' 'NF==5' "$j"_flank_analysis | awk -F '\t' '{print $0,$3"_"$5}' | awk '!seen[$6]++' | cat -n | sed "s/^ *//g" | sed "s/\t/,/g" | cut -f1 -d " ")
do
linenumber=$(echo $i | cut -f1 -d ',')
genome_name=$(echo $i | cut -f2 -d "@" | cut -f1 -d "_")
ortholog_name=$(echo $i | cut -f2 -d ',')
grep -w $ortholog_name ../nr_gtf/$genome_name.filtered.nr.gtf | sed "s/ \"/ \""$j"_"$linenumber"_/g" >> proxflanks_interim.gtf #Get a concatenated gtf
done
done

#Get the sequences
#The $4>0 is part is to get rid of those cases in which the start codon dips to negative
#I'm making the choice to skip those that don't have 500bp on either edge
cat proxflanks_interim.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-500,$4,$6,$7,$8,$9}' | awk -F '\t' '($4>0)' | sed "s/ \"/ \"left500_/g" | gtf2bed | bedtools getfasta -s -name -fi ../Ecoli.genome.pangenomedb.fna -bed - > proxflanks.fna
cat proxflanks_interim.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+500,$6,$7,$8,$9}' | sed "s/ \"/ \"right500_/g" | gtf2bed | bedtools getfasta -s -name -fi ../Ecoli.genome.pangenomedb.fna -bed - >> proxflanks.fna

#Remove those that lack either or both flanks:
#The code below counts the number of commas per line after massaging it, and if there are no commas——that means at least one flank is missing
#Is the following even needed? Checks if we have 500bp on both ends...
#grep "^>" proxflanks.fna | cut -f1 -d ":" | sed "s/_/\t/" | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | grep "," | cut -f1 | grep --no-group-separator -A1 -f - proxflanks.fna > interim
#mv interim proxflanks.faa

#First, the blast searches (Time consuming)
blastn -query geneflanks.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 109734 -max_hsps 1 -out geneflanks_extragenus_blastn
#109734
blastn -query geneflanks.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Ecoli_excluded.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 195119 -max_hsps 1 -out geneflanks_intragenus_blastn
#195119
blastn -query geneflanks.faa -db ../Ecoli.genome.pangenomedb -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 1761 -max_hsps 1 -out geneflanks_pangenome_blastn
#1761
blastn -query proxflanks.fna -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 109734 -max_hsps 1 -out proxflanks_extragenus_blastn
blastn -query proxflanks.fna -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Ecoli_excluded.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 195119 -max_hsps 1 -out proxflanks_intragenus_blastn
blastn -query proxflanks.fna -db ../Ecoli.genome.pangenomedb -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 1761 -max_hsps 1 -out proxflanks_pangenome_blastn
#choice: max-hsps 1

#Do next steps after blasts are done
