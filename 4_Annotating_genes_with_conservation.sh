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

#Goal:

1. Assign each cluster a value of taxonomic restriction
2. Assign each gene a cluster
3. Each gff should have a value assigned to it

#1 is toughest, possibly
#First, identify the contigs in which they have hits
cat all_proteins_vs_GBRS_annotated.tsv_xaa all_proteins_vs_GBRS_annotated.tsv_xab all_proteins_vs_GBRS_annotated.tsv_xac all_proteins_vs_GBRS_annotated.tsv_xad | awk -F '\t' '($5>60&&$16<0.001)' | cut -f1,2 > GBRS_annotated_hits.interim.tsv
awk -F '\t' '($5>60&&$16<0.001)' all_proteins_vs_GBRS_ORFs.tsv | cut -f1,2 > GBRS_ORF_hits.interim.tsv

############
############
############
############

#Pangenome analysis:

cut -f1 -d "@" 041124_all_protein_clusters.tsv | sed "s/\t/\tStrain_/g" | awk '
{
    row[$1];           # Store unique row (column 1) values
    col[$2];           # Store unique column (column 2) values
    matrix[$1][$2]=1;  # Mark presence of value in the matrix
}
END {
    # Print header row
    printf "    ";
    for (c in col) printf c "\t"; 
    print "";

    # Print each row, with 0s and 1s
    for (r in row) {
        printf r "\t";
        for (c in col) {
            if (matrix[r][c]) printf "1\t"; 
            else printf "0\t";
        }
        print "";
    }
}' > 041124_all_protein_clusters.wide.tsv

sed "s/\t/,/g" 041124_all_protein_clusters.wide.tsv > 041124_all_protein_clusters.wide.csv
#nano to add in gene in first line
python3 parsing2.py 041124_all_protein_clusters.wide.csv > 041124_all_protein_clusters.wide.mod.csv

rev 041124_all_protein_clusters.wide.mod.csv | cut -f2- -d ',' | rev > test && mv test 041124_all_protein_clusters.wide.mod.csv

#First make files containing just singletons and strict cores (which are easily extractable)
awk -F'Strain' 'NF>1' 041124_all_protein_clusters.wide.mod.csv | #removes all genes that are absent from all 469
tail -n+3 | #Remove header and "lineage" line
awk -F'Strain' 'NF==2' | 
cut -f1 -d ',' | sort -u > 419_singletons.txt

#Strict cores:

grep -v ",0" 041124_all_protein_clusters.wide.mod.csv |
tail -n+2 | #The "lineage" line contains ,0s, so auto removed in the last step
cut -f1 -d ',' | sort -u > 419_strictcores.txt

grep -wf ../419_isolate_names.txt ../500_ipp_lineagedesignations.tsv > 419_isolates_lineages.tsv

#Per the scheme below, phylogroups are identified as "@integer@", and lineages are "@integer%"
awk '{OFS="\t"}{print "Strain_"$1,"Strain_"$1"@"$4"@"$5"%"}' 419_isolates_lineages.tsv > 419_strain_lineage_classifications.tsv

#The phylogroup and lineage classifications for each strain are then added on to the p/a table:
python3 parsing3.py 419_strain_lineage_classifications.tsv 041124_all_protein_clusters.wide.mod.csv > 041124_all_protein_clusters.lineage.presenceabsence.csv

#We then take out strict cores, singletons, complete blanks from the resultant file:

awk -F'Strain' 'NF>2' 041124_all_protein_clusters.lineage.presenceabsence.csv | #singletons and blanks removed
grep ",0" |  #At least one strain that lacks the gene - i.e., strictcores removed
tail -n+2 > 041124_all_protein_clusters.lineage.presenceabsence.noncore.nonsingle.csv

#We now convert this file to a matrix, reporting #strains per phylogroup/lineage containing each gene
#The two lines of code below achieve this for phylogroups and lineages, respectively:

python3 count_strings_phylogroups.py 041124_all_protein_clusters.lineage.presenceabsence.noncore.nonsingle.csv > 041124_all_protein_clusters.PresenceAbsencePhylogroups.csv
python3 count_strings_lineages.py 041124_all_protein_clusters.lineage.presenceabsence.noncore.nonsingle.csv > 041124_all_protein_clusters.PresenceAbsenceLineages.csv

#Convert raw counts to a percentage of #strain per group

cut -f4 419_isolates_lineages.tsv | sort | uniq -c | sed "s/^ *//g" | sed "s/ /,/g" | awk -F ',' '{OFS=FS}{print "@"$2"@",$1}'> phylogroup_totalstrains.csv
cut -f5 419_isolates_lineages.tsv | sort | uniq -c | sed "s/^ *//g" | sed "s/ /,/g" | awk -F ',' '{OFS=FS}{print "@"$2"%",$1}'> lineage_totalstrains.csv

#Transform raw p/a values to those normalized by (= divided by) total strains per clade

python3 normalize_gene_counts.py 041124_all_protein_clusters.PresenceAbsencePhylogroups.csv phylogroup_totalstrains.csv | sed '/^$/d' > PresenceAbsencePhylogroups.transformed.csv
python3 normalize_gene_counts_lineages.py 041124_all_protein_clusters.PresenceAbsenceLineages.csv lineage_totalstrains.csv | sed '/^$/d' > PresenceAbsenceLineages.transformed.csv

#Phylogroup-specific core genes
python3 phylogroup_gain.py PresenceAbsencePhylogroups.transformed.csv > phylogroup_gain.csv
python3 phylogroup_gain.py PresenceAbsenceLineages.transformed.csv > lineage_gain.csv

#Near-core cases:

awk -F',0' 'NF==2' 041124_all_protein_clusters.lineage.presenceabsence.noncore.nonsingle.csv | cut -f1 -d ',' | sort -u > nearcore_1strainloss.txt
awk -F',0' 'NF==3' 041124_all_protein_clusters.lineage.presenceabsence.noncore.nonsingle.csv | cut -f1 -d ',' | sort -u > nearcore_2strainloss.txt
awk -F',0' 'NF==4' 041124_all_protein_clusters.lineage.presenceabsence.noncore.nonsingle.csv | cut -f1 -d ',' | sort -u > nearcore_3strainloss.txt

#1-phylogroup loss:
python3 phylogroup_loss.py PresenceAbsencePhylogroups.transformed.csv > phylogroup_loss.csv
#1-Lineage loss:
python3 phylogroup_loss.py PresenceAbsenceLineages.transformed.csv > lineage_loss.csv
