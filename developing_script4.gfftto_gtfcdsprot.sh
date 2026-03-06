#Converting gffs to gtf:
#Prodigal:
for i in $(ls prodigal/*_prodigal.gff | cut -f1 -d "_" | cut -f2 -d "/")
do
grep -v "#" prodigal/"$i"_prodigal.gff |
awk -F '\t' '($3=="CDS")' |
cut -f-8 | cat -n | sed "s/^ *//g" |
awk -F '\t' '{OFS=FS}{print $2,$3,$4,$5,$6,$7,$8,$9,"transcript_id \"prod@"$2"@"$1"\";gene_id \"prod@"$2"@"$1"\";"}' > prodigal/"$i"_prodigal.gtf
done

#Genemarks2:
for i in $(ls genemarks2/*_genemarks2.gff | cut -f1 -d "_" | cut -f2 -d "/")
do
grep -v "#" genemarks2/"$i"_genemarks2.gff |
awk -F '\t' '($3=="CDS")' |
cut -f-8 | cat -n | sed "s/^ *//g" |
awk -F '\t' '{OFS=FS}{print $2,$3,$4,$5,$6,$7,$8,$9,"transcript_id \"gms2@"$2"@"$1"\";gene_id \"gms2@"$2"@"$1"\";"}' > genemarks2/"$i"_genemarks2.gtf
done

#gtf to cds:
#prodigal:
for i in $(ls prodigal/*gtf | cut -f1 -d "_" | cut -f2 -d "/")
do
cat prodigal/"$i"_prodigal.gtf | gtf2bed | bedtools getfasta -s -name -fi ../Ecoli_nonred_genomes/"$i"_genomic.renamed.fna -bed - > "$i"_prodigal.fna
done

#genemarks2:
for i in $(ls genemarks2/*gtf | cut -f1 -d "_" | cut -f2 -d "/")
do
cat genemarks2/"$i"_genemarks2.gtf | gtf2bed | bedtools getfasta -s -name -fi ../Ecoli_nonred_genomes/"$i"_genomic.renamed.fna -bed - > "$i"_genemarks2.fna
done

sed -i 's/::.*//' *.fna

#Make proteins:

for i in $(ls *fna | rev | cut -f2- -d "." | rev)
do
/stor/work/Ochman/hassan/tools/faTrans -stop "$i".fna "$i".prot.faa
done

#linearize protein multifastas:

for i in *prot.faa
do
seqkit fx2tab "$i" | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" > test && mv test "$i"
done

mkdir -p prodigal_cds
mkdir -p prodigal_prot
mkdir -p genemarks2_cds
mkdir -p genemarks2_prot

mv *prodigal.fna prodigal_cds
mv *genemarks2.fna genemarks2_cds
mv *prodigal.prot.faa prodigal_prot
mv *genemarks2.prot.faa genemarks2_prot

