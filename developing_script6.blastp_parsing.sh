#Outside genus, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_queries.prot.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.proteins.dmnd --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_GBRS_annotated.tsv -k 0 -b8 -c1
#Outside genus, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_queries.prot.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.genomes.orfipy.prot.dmnd --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_GBRS_ORFs.tsv -k 0 -b8 -c1
#Outside species, annotated 
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_queries.prot.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Ecoli_excluded.proteins.dmnd --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_annotated.tsv -k 0 -b8 -c1
#Outside species, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_queries.prot.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Ecoli_excluded.genomes.orfipy.prot.dmnd --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_ORFs.tsv -k 0 -b8 -c1
#Pangenome, annotated 
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_queries.prot.faa -d Ecoli.prot.pangenomedb.dmnd --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_annotated.tsv -k 0 -b8 -c1
#Pangenome, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_queries.prot.faa -d Ecoli.genome.pangenomedb.orfipy.dmnd --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_ORFs.tsv -k 0 -b8 -c1

cat all_proteins_vs_GBRS_annotated.tsv |
grep -v -F -f /stor/work/Ochman/hassan/MS_Ecoli_ORFans_Ch3/rethinking_clustering/all_Escherichia_slippedby_IDs.txt - | #remove contaminants
awk -F '\t' '($5>60&&$16<0.001)' | cut -f1 | sort -u > step1_genusspecific_nonORFan.txt

/stor/work/Ochman/hassan/MS_Ecoli_ORFans_Ch3/rethinking_clustering/all_Escherichia_slippedby_IDs.txt

#NOW TIME TO make presence/absence tables

#First we parse the blast files:

#Annotated, protein IDs:
awk -F '\t' '($16<0.001){print $1"\t"$2"\t"($5>60?"present":"partial")}' all_proteins_vs_GBRS_annotated.tsv \
| sort -u > gbrs.prot.pa.tsv

awk -F '\t' '($16<0.001){print $1"\t"$2"\t"($5>60?"present":"partial")}' all_proteins_vs_noncoliEscherichia_annotated.tsv \
| sed "s/(+)//g" | sed "s/(-)//g" | sort -u > intra.prot.pa.tsv

awk -F '\t' '($16<0.001){print $1"\t"$2"\t"($5>60?"present":"partial")}' all_proteins_vs_pangenome_annotated.tsv \
| sort -u > pan.prot.pa.tsv

# ORFs -> contig IDs
awk -F '\t' '($16<0.001){print $1"\t"$2"\t"($5>60?"present":"partial")}' all_proteins_vs_GBRS_ORFs.tsv \
| sed "s/_ORF/\t/g" | cut -f1,2,4 \
| sort -u > gbrs.contig.pa.tsv

awk -F '\t' '($16<0.001){print $1"\t"$2"\t"($5>60?"present":"partial")}' all_proteins_vs_noncoliEscherichia_ORFs.tsv \
| sed "s/_ORF/\t/g" | cut -f1,2,4 \
| sort -u > intra.contig.pa.tsv

awk -F '\t' '($16<0.001){print $1"\t"$2"\t"($5>60?"present":"partial")}' all_proteins_vs_pangenome_ORFs.tsv \
| sed "s/_ORF/\t/g" | cut -f1,2,4 \
| sort -u > pan.contig.pa.tsv

#We then join them with taxon info:

sort -k2 gbrs.prot.pa.tsv | join -1 2 -2 2 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_extragenus_protein_contig_taxa.reduced.noescherichia.tsv | awk '{OFS="\t"}{print $2,$4,$3,$5}' > diamond.presenceabsence.tsv1
sort -k2 intra.prot.pa.tsv | join -1 2 -2 2 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_intragenus_genome_protein_taxa.tsv | awk '{OFS="\t"}{print $2,$4,$3,$5}' >> diamond.presenceabsence.tsv2
sed "s/\t/@/g" pan.prot.pa.tsv | sed "s/_/@/2" | cut -f1,2,3,5,8 -d "@" | sed "s/@/\t/3" | sed "s/@/\t/3" | sort -k2 | join -1 2 -2 1 - ../genome_cluster.tsv | awk '{OFS="\t"}{print $2,$1,$3,$4}' > diamond.presenceabsence.tsv3
sort -k2 gbrs.contig.pa.tsv | join -1 2 -2 2 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_extragenus_genome_contig_taxa.reduced.noescherichia.tsv | awk '{OFS="\t"}{print $2,$4,$3,$5}' > diamond.presenceabsence.tsv4
sort -k2 intra.contig.pa.tsv | join -1 2 -2 2 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_intragenus_genome_contig_taxa.tsv | awk '{OFS="\t"}{print $2,$4,$3,$5}' > diamond.presenceabsence.tsv5
awk '{print $1,$3,$2}' pan.contig.pa.tsv | cut -f1,2 -d "_" | awk '{OFS="\t"}{print $1,$3,$2}' | sort -k2 | join -1 2 -2 1 - ../genome_cluster.tsv | awk '{OFS="\t"}{print $2,$1,$3,$4}' > diamond.presenceabsence.tsv6

#Collapse into one file:

cat *diamond.presenceabsence* | awk -F '\t' 'BEGIN{OFS="\t"}
{
    key = $1 SUBSEP $2 SUBSEP $4
    if (!(key in best) || $3 == "present")
        best[key] = $3
}
END {
    for (key in best) {
        split(key, a, SUBSEP)
        print a[1], a[2], best[key], a[3]
    }
}' | sort -t $'\t' -k1,1 -k2,2 -k4,4 > all_genes_of_interest.presence.tsv

