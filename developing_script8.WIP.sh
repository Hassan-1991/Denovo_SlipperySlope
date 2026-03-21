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
for j in $(cut -f1 -d "(" ../step1_genusspecific_ORFan.replaced.txt)
do
    echo "#!/bin/bash" > "${j}.sh"

    echo "for i in \$(awk '(\$1==\"$j\")' Ecoli.nr.filtered.prot.clusters.longestrepresentative.silixcollapsed.tsv | cut -f2)" >> "${j}.sh"
    echo "do" >> "${j}.sh"
    echo "    genome_name=\$(echo \$i | cut -f2 -d \"@\" | cut -f1 -d \"_\")" >> "${j}.sh"
    echo "    gene_ortholog_gtf=\$(grep -w \$i ../nr_gtf/\$genome_name.filtered.nr.gtf)" >> "${j}.sh"
    echo "    echo \$i > ${j}_temp" >> "${j}.sh"
    echo "    echo \$gene_ortholog_gtf | sed \"s/ /%/9\" | sed \"s/ /%/9\" | sed \"s/ /\\t/g\" | sed \"s/%/ /g\" | bedtools closest -a - -b ../nr_gtf/\$genome_name.filtered.nr.gtf -D a -id -io | awk -F '\\t' '(\$12==\"CDS\")' | cut -f 6 -d \"\\\"\" >> \"\$i\"_upstream" >> "${j}.sh"
    echo "    cat \"\$i\"_upstream >> ${j}_temp" >> "${j}.sh"
    echo "    grep -w -F -f \"\$i\"_upstream Ecoli.nr.filtered.prot.clusters.longestrepresentative.silixcollapsed.tsv | cut -f1 >> ${j}_temp" >> "${j}.sh"
    echo "    rm \"\$i\"_upstream" >> "${j}.sh"
    echo "    echo \$gene_ortholog_gtf | sed \"s/ /%/9\" | sed \"s/ /%/9\" | sed \"s/ /\\t/g\" | sed \"s/%/ /g\" | bedtools closest -a - -b ../nr_gtf/\$genome_name.filtered.nr.gtf -D a -iu -io | awk -F '\\t' '(\$12==\"CDS\")' | cut -f 6 -d \"\\\"\" >> \"\$i\"_downstream" >> "${j}.sh"
    echo "    cat \"\$i\"_downstream >> ${j}_temp" >> "${j}.sh"
    echo "    grep -w -F -f \"\$i\"_downstream Ecoli.nr.filtered.prot.clusters.longestrepresentative.silixcollapsed.tsv | cut -f1 >> ${j}_temp" >> "${j}.sh"
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
