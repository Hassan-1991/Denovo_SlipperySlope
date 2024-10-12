#Outside genus, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_proteins/04072024_nonEscherichia_protein --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_GBRS_annotated.tsv -k 0 -b8 -c1
#Takes a long time:
xaa - done
xab - done
xac - done
xad - ongoing (pod1)
#Outside genus, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/04072024_nonEscherichia_genomes.getorf.noCTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_GBRS_ORFs.tsv -k 0 -b8 -c1
done
#Outside species, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_annotated.tsv -k 0 -b8 -c1
done
#Outside species, ORFs_ATG
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes.getorf.ATG --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_ORFs_ATG.tsv -k 0 -b8 -c1
ongoing (pod2)
#Outside species, ORFs_TTG_GTG
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q all_proteins_reps.faa -d /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/10082024_noncoli_Escherichia_database/noncoliEscherichia_allgenomes.getorf.TTG_GTG --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_noncoliEscherichia_ORFs_TTG_GTG.tsv -k 0 -b8 -c1

#Salmonella ORF:
#ongoing:
grep --no-group-separator -A1 -f non-Salmonella_identifiers.txt /stor/scratch/Ochman/hassan/RefSeq_Complete_Genomes_04062024/04062024_RS_GB_complete_bacterial_genomes/all_GBRS_getorf.noCTG > 04072024_nonSalmonella_genomes.getorf.noCTG
#Still to run:
/stor/work/Ochman/hassan/tools/faTrans -stop 04072024_nonSalmonella_genomes.getorf.noCTG 04072024_nonSalmonella_genomes.getorf.noCTG.prot
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in 04072024_nonSalmonella_genomes.getorf.noCTG.prot --db 04072024_nonSalmonella_genomes.getorf.noCTG.prot
