#Concatenate, proteins:
cat prodigal_prot/*prot.faa genemarks2_prot/*prot.faa > Ecoli.prot.pangenomedb.faa

#Concatenate, genomes:
cat ../Ecoli_nonred_genomes/*_genomic.renamed.fna > Ecoli.genome.pangenomedb.fna

#Extract proteins encoded from all ORFs:
conda activate orfipy
orfipy --start ATG,TTG,GTG --stop TGA,TAG,TAA --pep Ecoli.genome.pangenomedb.orfipy Ecoli.genome.pangenomedb.fna

#Arrange everything in directories

#Make protein and genome databases
makeblastdb -in Ecoli.genome.pangenomedb.fna -dbtype nucl -out Ecoli.genome.pangenomedb
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Ecoli.prot.pangenomedb.faa --db Ecoli.prot.pangenomedb
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in orfipy_Ecoli.genome.pangenomedb.fna_out/Ecoli.genome.pangenomedb.orfipy --db Ecoli.genome.pangenomedb.orfipy

#Get back to the part above later

#Non-redundantify and filter genomes

#Genes located at edges:

for d in prodigal genemarks2; do
  for f in "$d"/*.gtf; do
    awk -F'\t' '
      !seen[$1]++ {first[$1]=$0; order[++n]=$1}
      {last[$1]=$0}
      END {
        for (j=1; j<=n; j++) {
          k=order[j]
          split(first[k], a, "\""); print a[2]
          if (first[k] != last[k]) {
            split(last[k], b, "\""); print b[2]
          }
        }
      }
    ' "$f"
  done > "${d}_edge_genes.txt"
done

#Lacking dual annotation:
#First, merge the two annotations (retain longer ORF):
mkdir -p nr_gtf

for f in prodigal/*_prodigal.gtf; do
  i=$(basename "$f" _prodigal.gtf)
  awk '{
    stop = ($7 == "+") ? $5 : $4
    len  = $5 - $4
    key  = $1 ":" $7 ":" stop
    if (!(key in best) || len > best[key]) {
      best[key] = len
      line[key] = $0
    }
  }
  END {
    for (k in line) print line[k]
  }' "$f" "genemarks2/${i}_genemarks2.gtf" > "nr_gtf/${i}.nr.gtf"
done

#Genes with single annotation:
: > annotcomp.tsv

for f in nr_gtf/*.nr.gtf; do
  i=$(basename "$f" .nr.gtf)

  awk -F'\t' '
    function get_gene_id(attr,   m) {
      if (match(attr, /gene_id "[^"]+"/)) {
        m = substr(attr, RSTART, RLENGTH)
        sub(/^gene_id "/, "", m)
        sub(/"$/, "", m)
        return m
      }
      return ""
    }

    FNR==NR {
      stop = ($7=="+") ? $5 : $4
      id = get_gene_id($9)
      key = $1 FS $7 FS stop
      nr[key] = id
      next
    }

    {
      stop = ($7=="+") ? $5 : $4
      key = $1 FS $7 FS stop
      if (key in nr) {
        id = nr[key]
        src = (FILENAME ~ /prodigal/) ? "prodigal" : "genemarks2"
        hits[id,src] = 1
      }
    }

    END {
      for (k in nr) {
        id = nr[k]
        out=""
        if (hits[id,"prodigal"])   out=out"prodigal,"
        if (hits[id,"genemarks2"]) out=out"genemarks2,"
        sub(/,$/,"",out)
        if (out!="") print id "\t" out
      }
    }
  ' "$f" "prodigal/${i}_prodigal.gtf" "genemarks2/${i}_genemarks2.gtf" >> annotcomp.tsv
done

awk -F '\t' '($2!~",")' annotcomp.tsv  | cut -f1 | sort -u > genes_with_oneannotation.txt

#Genes that begin with other than ATG, TTG or GTG:

for d in prodigal genemarks2; do
  : > "${d}_genes_without_start.txt"
  : > "${d}_genes_without_stop.txt"

  for f in "$d"/*.gtf; do
    i=$(basename "$f")
    i=${i%_${d}.gtf}

    bedtools getfasta -s -name \
      -fi "../Ecoli_nonred_genomes/${i}_genomic.renamed.fna" \
      -bed <(gtf2bed < "$f") |
    seqkit fx2tab |
    awk -F'\t' '
      {
        gene=$1
        seq=toupper($2)
        sub(/\(.*/, "", gene)

        start=substr(seq,1,3)
        stop=substr(seq,length(seq)-2,3)

        if (start!="ATG" && start!="TTG" && start!="GTG")
          print gene >> start_out

        if (stop!="TGA" && stop!="TAA" && stop!="TAG")
          print gene >> stop_out
      }
    ' start_out="${d}_genes_without_start.txt" stop_out="${d}_genes_without_stop.txt"

  done
done

sed -i 's/::.*//' *without*txt

#Now to exclude these from gtfs:

for d in prodigal genemarks2; do
  prefix=$([ "$d" = prodigal ] && echo '^prod' || echo '^gms2')

  cat "${d}_genes_without_stop.txt" \
      "${d}_genes_without_start.txt" \
      "${d}_edge_genes.txt" \
      genes_with_oneannotation.txt |
  grep "$prefix" | sort -u > "${d}_exclude.txt"

  for f in "$d"/*.gtf; do
    i=$(basename "$f" "_${d}.gtf")
    awk -F'\t' '
      NR==FNR {exclude[$1]=1; next}
      {
        if (match($9, /gene_id "[^"]+"/)) {
          id = substr($9, RSTART+9, RLENGTH-10)
          if (!(id in exclude)) print
        }
      }
    ' "${d}_exclude.txt" "$f" > "${d}/${i}_${d}.filtered.gtf"
  done
done

for f in prodigal/*_prodigal.filtered.gtf; do
  i=$(basename "$f" _prodigal.filtered.gtf)
  awk '{
    stop = ($7 == "+") ? $5 : $4
    len  = $5 - $4
    key  = $1 ":" $7 ":" stop
    if (!(key in best) || len > best[key]) {
      best[key] = len
      line[key] = $0
    }
  }
  END {
    for (k in line) print line[k]
  }' "$f" "genemarks2/${i}_genemarks2.filtered.gtf" > "nr_gtf/${i}.filtered.nr.gtf"
done

####RUN FROM HERE TOMORROW:::#####

#Extract CDS:

for i in $(cat 460_genomes.txt)
do
cat Ecoli.460.gtfs/"$i".filtered.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi Ecoli.460.genomes/"$i" -bed - > Ecoli.460.CDS/"$i".filtered.CDS.faa
done

#Translate to protein:

for i in $(cat 460_genomes.txt)
do
/stor/work/Ochman/hassan/tools/faTrans -stop Ecoli.460.CDS/"$i".filtered.CDS.faa Ecoli.460.prots/"$i".filtered.prot.faa
done

#Concatenate:
for i in $(cat 460_genomes.txt)
do
genename=$(echo $i | cut -f1,2 -d "_")
sed "s/>/>"$genename"@/g" Ecoli.460.prots/"$i".filtered.prot.faa >> Ecoli.all.filtered.prot.faa
sed "s/>/>"$genename"@/g" Ecoli.460.CDS/"$i".filtered.CDS.faa >> Ecoli.all.filtered.CDS.faa
done

#Number of proteins - 2,179,313
mkdir mmseqs_covmode0
cp Ecoli.all.filtered.prot.faa mmseqs_covmode0
cd mmseqs_covmode0

mmseqs createdb Ecoli.all.filtered.prot.faa Ecoli.all.filtered.prot

mmseqs search Ecoli.all.filtered.prot Ecoli.all.filtered.prot resultDB tmp --min-seq-id 0.8 -c 0.8 --cov-mode 0
mmseqs convertalis Ecoli.all.filtered.prot Ecoli.all.filtered.prot resultDB resultDB.m8
mmseqs linclust Ecoli.all.filtered.prot clusterDB tmp --min-seq-id 0.8 -c 0.8 --cov-mode 0
mmseqs createtsv Ecoli.all.filtered.prot Ecoli.all.filtered.prot clusterDB Ecoli.all.filtered.prot.clusters.covmode0.tsv

#Extract longest sequence from each cluster
#For this purpose, first extract lengths
seqkit fx2tab Ecoli.all.filtered.prot.faa | awk -F '\t' '{print $1,length($2)}' | sort -k1 > Ecoli.all.filtered.prot.lengths.tsv

sort -k2 mmseqs_covmode0/Ecoli.all.filtered.prot.clusters.covmode0.tsv | join -1 2 -2 1 - Ecoli.all.filtered.prot.lengths.tsv | #attaching lengths to each gene
sort -k2 | awk '{print $0 | "sort -k2,2 -k3,3n"}' | awk '!seen[$2]++' | sort -k2 > representative_sequences.interim.tsv #Longest sequence per cluster retained

#Make a three-column tsv file which lists longest sequence per cluster in col1, representative sequence picked by mmseqs2 in column3, and genes belonging to that cluster in col2
sort -k1 mmseqs_covmode0/Ecoli.all.filtered.prot.clusters.covmode0.tsv | join -1 1 -2 2 -o 1.1 2.2 1.2 2.1 - representative_sequences.interim.tsv | awk '{print $4,$3,$1}' | sed "s/ /\t/g" > Ecoli.all.filtered.prot.clusters.longestrepresentative.tsv

#Get CDS and protein sequences
seqkit fx2tab Ecoli.all.filtered.prot.faa | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" > Ecoli.all.filtered.prot.linear.faa
cut -f1 Ecoli.all.filtered.prot.clusters.longestrepresentative.tsv | sort -u | grep --no-group-separator -A1 -F -f - Ecoli.all.filtered.prot.linear.faa > Ecoli.all.filtered.prot.clusters.longestrepresentative.faa
cut -f1 Ecoli.all.filtered.prot.clusters.longestrepresentative.tsv | sort -u | grep --no-group-separator -A1 -F -f - Ecoli.all.filtered.CDS.faa > Ecoli.all.filtered.prot.clusters.longestrepresentative.CDS.faa

#Extract the sequences in each cluster and place them in files
cut -f1 mmseqs2_recluster_silix_output.tsv | sort | uniq -c | grep -v " 1 " | rev | cut -f1 -d " " | rev | sed "s/$/\t/g" | grep -F -f - mmseqs2_recluster_silix_output.tsv | sed "s/\t/,/g" > mmseqs2_recluster_silix_output.interim.tsv
for i in $(cat mmseqs2_recluster_silix_output.interim.tsv | cut -f1 -d ",")
do
echo $i | sed "s/$/,/g" | grep -F -f - mmseqs2_recluster_silix_output.interim.tsv | cut -f2 -d ',' | sed "s/^/>/g" | grep --no-group-separator -A1 -F -f - ../Ecoli.all.filtered.prot.clusters.longestrepresentative.faa > "$i".seqs.faa
done

#replace:
grep "^>" Ecoli.all.filtered.prot.clusters.longestrepresentative.faa | tr -d ">" > Ecoli.all.filtered.prot.clusters.longestrepresentative.txt
awk 'FNR==NR {rep[$2]=$1; next} {for(i=1;i<=NF;i++) if ($i in rep) $i=rep[$i]; print}' mmseqs2_silix/replacements.tsv Ecoli.all.filtered.prot.clusters.longestrepresentative.txt | sort -u > Ecoli.all.filtered.prot.clusters.longestrepresentative.replaced.txt

grep --no-group-separator -A1 -F -f Ecoli.all.filtered.prot.clusters.longestrepresentative.replaced.txt Ecoli.all.filtered.prot.linear.faa > Ecoli.all.filtered.prot.clusters.longestrepresentative.replaced.faa
