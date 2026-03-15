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
