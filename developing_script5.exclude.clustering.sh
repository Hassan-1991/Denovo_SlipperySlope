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

mkdir -p nr_CDS
mkdir -p nr_prot

#Extract CDS:

for i in $(ls nr_gtf/*filtered.nr.gtf | rev | cut -f1 -d "/" | rev | cut -f1 -d ".")
do
cat nr_gtf/"$i".filtered.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi ../Ecoli_nonred_genomes/"$i"_genomic.renamed.fna -bed - > "$i".filtered.nr.CDS.fna
done

sed -i 's/::.*//' *filtered.nr.CDS.fna

#Translate to protein:

for i in $(ls nr_gtf/*filtered.nr.gtf | rev | cut -f1 -d "/" | rev | cut -f1 -d ".")
do
/stor/work/Ochman/hassan/tools/faTrans -stop "$i".filtered.nr.CDS.fna "$i".filtered.nr.prot.faa
done

for i in *filtered.nr.prot.faa
do
seqkit fx2tab "$i" | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" > test && mv test "$i"
done

mv *filtered.nr.CDS.fna nr_CDS
mv *filtered.nr.prot.faa nr_prot

#Concatenate:
cat nr_CDS/* > Ecoli.nr.filtered.CDS.faa
cat nr_prot/* > Ecoli.nr.filtered.prot.faa

#Number of proteins - 2,420,394
mkdir -p mmseqs_covmode0
cp Ecoli.nr.filtered.prot.faa mmseqs_covmode0
cd mmseqs_covmode0

mmseqs createdb Ecoli.nr.filtered.prot.faa Ecoli.nr.filtered.prot
mmseqs linclust Ecoli.nr.filtered.prot clusterDB tmp --min-seq-id 0.6 -c 0.6 --cov-mode 0
mmseqs createtsv Ecoli.nr.filtered.prot Ecoli.nr.filtered.prot clusterDB Ecoli.nr.filtered.prot.clusters.covmode0.tsv

#Extract longest sequence from each cluster
seqkit fx2tab Ecoli.nr.filtered.prot.faa | awk -F '\t' '{print $1"\t"length($2)}' | sort -k1 > Ecoli.nr.filtered.prot.lengths.tsv

# make longest-member / member / mmseqs-representative table
awk '
BEGIN{FS=OFS="\t"}
NR==FNR {
    len[$1]=$2
    next
}
{
    rep=$1
    mem=$2

    if (!(rep in longest_len) || len[mem] > longest_len[rep]) {
        longest_len[rep]=len[mem]
        longest_id[rep]=mem
    }

    reps[NR]=rep
    mems[NR]=mem
    n=NR
}
END{
    for (i=1; i<=n; i++) {
        print longest_id[reps[i]], mems[i], reps[i]
    }
}
' Ecoli.nr.filtered.prot.lengths.tsv Ecoli.nr.filtered.prot.clusters.covmode0.tsv \
> Ecoli.nr.filtered.prot.clusters.longestrepresentative.tsv

#Get CDS and protein sequences
#clustering done, so exit directory
cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_genomics/Ecoli.clusters.0.1.gaps/Ecoli_nonred_gffs
cp mmseqs_covmode0/Ecoli.nr.filtered.prot.clusters.longestrepresentative.tsv .
cut -f1 Ecoli.nr.filtered.prot.clusters.longestrepresentative.tsv | sort -u | sed "s/^/>/g" | grep --no-group-separator -w -A1 -F -f - Ecoli.nr.filtered.prot.faa > Ecoli.nr.filtered.prot.longestrepresentative.prot.faa
cut -f1 Ecoli.nr.filtered.prot.clusters.longestrepresentative.tsv | sort -u | sed "s/^/>/g" | grep --no-group-separator -w -A1 -F -f - Ecoli.nr.filtered.CDS.faa > Ecoli.nr.filtered.prot.longestrepresentative.CDS.fna

###Silix:

/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Ecoli.nr.filtered.prot.longestrepresentative.prot.faa --db Ecoli.nr.filtered.prot.longestrepresentative.prot
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli.nr.filtered.prot.longestrepresentative.prot.faa -d Ecoli.nr.filtered.prot.longestrepresentative.prot --outfmt 6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore --ultra-sensitive --out Ecoli.nr.filtered.prot.longestrepresentative.prot.allvall.tsv -k 0 -b8 -c1
awk -F '\t' '($3>60&&$4>60&&$12<0.001)' Ecoli.nr.filtered.prot.longestrepresentative.prot.allvall.tsv | cut -f1-3,5- > silix_input.tsv
/stor/work/Ochman/hassan/tools/silix-1.3.0/src/silix -f cluster_ -i 0.6 -r 0.6 -q -2 -s 3 Ecoli.nr.filtered.prot.longestrepresentative.prot.faa silix_input.tsv | sort -k2 > silix_output.tsv

#Replace the mmseqs2 cluster representatives with the silix cluster reps:

awk '
BEGIN{
    FS=OFS="\t"   # input and output are tab-delimited
}

# First file: lengths table
# Store length of each protein ID
NR==FNR{
    len[$1] = $2
    next
}

# Second file: silix_output.tsv
# col1 = silix cluster ID
# col2 = protein ID
{
    silix_cluster = $1
    prot          = $2

    # If this is the first protein seen for this SiLiX cluster,
    # or this protein is longer than the current best one,
    # then update the "best representative" for that SiLiX cluster
    if (!(silix_cluster in best_len) || len[prot] > best_len[silix_cluster]) {
        best_len[silix_cluster] = len[prot]
        best_rep[silix_cluster] = prot
    }

    # Also remember which SiLiX cluster each protein belongs to
    prot2cluster[prot] = silix_cluster
}

END{
    # Output:
    # col1 = protein ID from silix_output.tsv
    # col2 = SiLiX cluster ID
    # col3 = longest protein in that SiLiX cluster
    for (prot in prot2cluster) {
        clus = prot2cluster[prot]
        print prot, clus, best_rep[clus]
    }
}
' mmseqs_covmode0/Ecoli.nr.filtered.prot.lengths.tsv silix_output.tsv \
| sort -k1,1 > silix_longestrep_map.tsv

#Replace

awk '
BEGIN{
    FS=OFS="\t"   # tab-delimited input/output
}

# First file: silix_longestrep_map.tsv
# col1 = old representative
# col2 = silix cluster
# col3 = new longest representative for that silix cluster
NR==FNR{
    replacement[$1] = $3
    next
}

# Second file: Ecoli.nr.filtered.prot.clusters.longestrepresentative.tsv
# col1 = current representative
# col2 = member
# col3 = mmseqs representative
{
    old_rep = $1

    # If this representative belongs to a SiLiX cluster,
    # replace it with the longest representative of that SiLiX cluster
    if (old_rep in replacement) {
        $1 = replacement[old_rep]
    }

    # Print updated row
    print
}
' silix_longestrep_map.tsv Ecoli.nr.filtered.prot.clusters.longestrepresentative.tsv \
> Ecoli.nr.filtered.prot.clusters.longestrepresentative.silixcollapsed.tsv

#Get prot and CDS:
cut -f1 Ecoli.nr.filtered.prot.clusters.longestrepresentative.silixcollapsed.tsv | sort -u | grep -w --no-group-separator -A1 -F -f - Ecoli.nr.filtered.prot.longestrepresentative.prot.faa > Ecoli_queries.prot.faa
cut -f1 Ecoli.nr.filtered.prot.clusters.longestrepresentative.silixcollapsed.tsv | sort -u | grep -w --no-group-separator -A1 -F -f - Ecoli.nr.filtered.prot.longestrepresentative.CDS.fna > Ecoli_queries.CDS.fna




