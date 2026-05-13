#These are steps that need to be incorporated way earlier in the pipeline
#Reconsider later

#Iterative clustering:

/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli.nr.filtered.prot.longestrepresentative.prot.faa -d Ecoli.nr.filtered.prot.longestrepresentative.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore --ultra-sensitive --out Ecoli.nr.filtered.prot.longestrepresentative.silixinput.tsv -k 0 -b8 -c1
/stor/work/Ochman/hassan/tools/silix-1.3.0/src/silix -f cluster_ -i 0.6 -r 0.6 -q -2 -s 3 Ecoli.nr.filtered.prot.longestrepresentative.prot.faa Ecoli.nr.filtered.prot.longestrepresentative.silixinput.tsv | sort -k2 > silix_output.tsv
seqkit fx2tab Ecoli.nr.filtered.prot.longestrepresentative.prot.faa | awk '{print $1,length($2)}' | sort -k1 | join -1 1 -2 2 - silix_output.tsv | sort -k3,3 -k2,2nr | awk '!seen[$3]++' | cut -f1 -d " " | grep --no-group-separator -w -A1 -F -f - Ecoli.nr.filtered.prot.longestrepresentative.prot.faa > Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa

#qcov1:

/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa --db Ecoli.nr.filtered.prot.longestrepresentative.reduced1
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa -d Ecoli.nr.filtered.prot.longestrepresentative.reduced1 --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore --ultra-sensitive --out Ecoli.nr.filtered.prot.longestrepresentative.reduced1.silixinput.tsv -k 0 -b8 -c1
/stor/work/Ochman/hassan/tools/silix-1.3.0/src/silix -f cluster_ -i 0.6 -r 0.6 -q 1 -s 3 Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa Ecoli.nr.filtered.prot.longestrepresentative.reduced1.silixinput.tsv | sort -k2 > silix_output.reduced1.tsv
seqkit fx2tab Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa | awk '{print $1,length($2)}' | sort -k1 | join -1 1 -2 2 - silix_output.reduced1.tsv | sort -k3,3 -k2,2nr | awk '!seen[$3]++' | cut -f1 -d " " | grep --no-group-separator -w -A1 -F -f - Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa > Ecoli.nr.filtered.prot.longestrepresentative.reduced2.faa

#qcov2

/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Ecoli.nr.filtered.prot.longestrepresentative.reduced2.faa --db Ecoli.nr.filtered.prot.longestrepresentative.reduced2
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli.nr.filtered.prot.longestrepresentative.reduced2.faa -d Ecoli.nr.filtered.prot.longestrepresentative.reduced2 --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore --ultra-sensitive --out Ecoli.nr.filtered.prot.longestrepresentative.reduced2.silixinput.tsv -k 0 -b8 -c1
/stor/work/Ochman/hassan/tools/silix-1.3.0/src/silix -f cluster_ -i 0.6 -r 0.6 -q 1 -s 3 Ecoli.nr.filtered.prot.longestrepresentative.reduced2.faa Ecoli.nr.filtered.prot.longestrepresentative.reduced2.silixinput.tsv | sort -k2 > silix_output.reduced2.tsv
seqkit fx2tab Ecoli.nr.filtered.prot.longestrepresentative.reduced2.faa | awk '{print $1,length($2)}' | sort -k1 | join -1 1 -2 2 - silix_output.reduced2.tsv | sort -k3,3 -k2,2nr | awk '!seen[$3]++' | cut -f1 -d " " | grep --no-group-separator -w -A1 -F -f - Ecoli.nr.filtered.prot.longestrepresentative.reduced2.faa > Ecoli.nr.filtered.prot.longestrepresentative.reduced3.faa

#Two different reduced protein candidate sets:

1. Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa
2. Ecoli.nr.filtered.prot.longestrepresentative.reduced2.faa

#Search against nr:

/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp \
  -q Ecoli.nr.filtered.prot.longestrepresentative.reduced1.faa \
  -d /stor/scratch/Ochman/hassan/nr/03272024_nr.dmnd \
  --taxonmap /stor/scratch/Ochman/hassan/nr/prot.accession2taxid.FULL \
  --taxonnames /stor/scratch/Ochman/hassan/nr/names.dmp \
  --taxonnodes /stor/scratch/Ochman/hassan/nr/nodes.dmp \
  --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore staxids sscinames \
  --ultra-sensitive \
  --out Ecoli.nr.filtered.prot.longestrepresentative.reduced1_vs_nr_withtax.tsv \
  -k 0 -b8 -c1

#Still ongoing

#parse:
cp /stor/scratch/Ochman/hassan/nr/nodes.dmp .

#Python code:

#!/usr/bin/env python3
import sys
from collections import defaultdict, deque

parent_to_children = defaultdict(list)

with open("nodes.dmp") as f:
    for line in f:
        parts = [x.strip() for x in line.split("|")]
        taxid, parent = parts[0], parts[1]
        parent_to_children[parent].append(taxid)

root = sys.argv[1]   # e.g. 562

q = deque([root])
while q:
    t = q.popleft()
    print(t)
    q.extend(parent_to_children[t])

#######

#!/usr/bin/env python3
import sys

nodes = "nodes.dmp"
taxid = sys.argv[1]

parent = {}

with open(nodes) as f:
    for line in f:
        parts = [x.strip() for x in line.split("|")]
        child, par = parts[0], parts[1]
        parent[child] = par

t = taxid
while True:
    print(t)
    if t == parent[t]:   # root
        break
    t = parent[t]

############

python get_ancestor_taxids.py 562 > Ecoli_562_ancestor_taxids.txt
python get_ancestor_taxids.py 620 >> Ecoli_562_ancestor_taxids.txt
sort -u Ecoli_562_ancestor_taxids.txt -o Ecoli_562_ancestor_taxids.txt
python get_descendant_taxids.py 562 > Ecoli_562_descendant_taxids.txt
python get_descendant_taxids.py 620 >> Ecoli_562_descendant_taxids.txt
python get_descendant_taxids.py 10239 >> Ecoli_562_descendant_taxids.txt
sort -u Ecoli_562_descendant_taxids.txt -o Ecoli_562_descendant_taxids.txt

awk -F '\t' '
ARGIND==1 { anc[$1]; next }
ARGIND==2 { desc[$1]; next }

{
  lname = tolower($NF)

  # drop rows with missing or specific unwanted labels
  if ($NF == "N/A" || lname ~ /uncultured[[:space:]]+escherichia|escherichia[[:space:]]+sp\./) next

  taxfield = $(NF-1)

  # ancestors: strict whole-field match
  if (taxfield in anc) next

  # descendants: loose token match
  n = split(taxfield, a, /[^0-9]+/)
  for (i=1; i<=n; i++)
    if (a[i] in desc) next

  print
}
' Ecoli_562_ancestor_taxids.txt \
  Ecoli_562_descendant_taxids.txt \
  Ecoli.nr.filtered.prot.longestrepresentative.reduced1_vs_nr_withtax.tsv > Ecoli.nr.filtered.prot.longestrepresentative.reduced1_vs_nr_withtax.badids_removed.tsv

awk -F '\t' '
($16 < 0.001) {
    OFS = FS
    status = ($5 > 60 ? "present" : "partial")
    print $1, status, $(NF-1), $NF
}
' Ecoli.nr.filtered.prot.longestrepresentative.reduced1_vs_nr_withtax.badids_removed.tsv > Ecoli.nr.filtered.prot.longestrepresentative.vs_nr.tsv

#Identify different classes of non-ORFans from this list:

#!/bin/bash

in="Ecoli.nr.filtered.prot.longestrepresentative.vs_nr.tsv"

awk -F'\t' -v OFS='\t' '
{
    print $1 > "all_col1.tmp"

    key = $1
    for (i = 3; i <= NF; i++) key = key OFS $i
    print key > "all_uniquehits.tmp"

    if ($2 == "present") {
        print $1 > "present_col1.tmp"
        print key > "present_uniquehits.tmp"
    }
}
' "$in"

sort -u all_col1.tmp > anyhit_nr_nonORFans.partial.txt

sort -u all_uniquehits.tmp |
cut -f1 |
sort |
uniq -c |
awk '($1 > 1)' > twohit_nr_nonORFans.partial.txt

sort -u present_col1.tmp > anyhit_nr_nonORFans.present.txt

sort -u present_uniquehits.tmp |
cut -f1 |
sort |
uniq -c |
awk '($1 > 1)' > twohit_nr_nonORFans.present.txt

rm -f all_col1.tmp all_uniquehits.tmp present_col1.tmp present_uniquehits.tmp

###
###
###
###
###
###

awk -F '\t' '
ARGIND==1 { anc[$1]; next }
ARGIND==2 { desc[$1]; next }

{
  lname = tolower($NF)

  # drop rows with missing or specific unwanted labels
  if ($NF == "N/A" || lname ~ /uncultured[[:space:]]+escherichia|escherichia[[:space:]]+sp\./) next

  taxfield = $(NF-1)

  # ancestors: strict whole-field match
  if (taxfield in anc) next

  # descendants: loose token match
  n = split(taxfield, a, /[^0-9]+/)
  for (i=1; i<=n; i++)
    if (a[i] in desc) next

  print
}
' empty.txt \
  virus_descendant_taxids.txt \
  Ecoli.nr.filtered.prot.longestrepresentative.reduced1_vs_nr_withtax.badids_removed.tsv > Ecoli.nr.filtered.prot.longestrepresentative.reduced1_vs_nr_withtax.badids_removed.2.tsv

#Gotta work in how to remove viral ids too

cat Ecoli.nr.filtered.prot.longestrepresentative.reduced1_vs_nr_withtax.tsv |
awk -F '\t' '($5>60&&$16<0.001)' |
awk -F '\t' '($19!~"Escherichia coli")' |
awk -F '\t' '($19!~"Shigella")' |
awk -F '\t' '($19!~"Escherichia sp.")' |
awk -F '\t' '($19!="N/A")' |
awk -F '\t' '($19!~"Enterobacteriaceae")' |
awk -F '\t' '($19!="Escherichia")' |
awk -F '\t' '($19!="Caudoviricetes sp.")' |
awk -F '\t' '($19!="unclassified Escherichia")' |
awk -F '\t' '($19!="Bacteriophage sp.")' |
awk -F '\t' '($19!="uncultured bacterium")' |
awk -F '\t' '($19!="Klebsiella pneumoniae IS22")' |
awk -F '\t' '($19!="Enterobacterales")' |
awk -F '\t' '($19!="Escherichia;uncultured bacterium")' |
awk -F '\t' '($19!="Enterobacter sp. EC-NT1")' |
awk -F '\t' '($19!="Escherichia;Caudoviricetes sp.")' |
awk -F '\t' '($19!="Salmonella sp. S13")' |
awk -F '\t' '($19!="Escherichia;Salmonella sp. S13")' |
awk -F '\t' '($19!~"phage")' |
cut -f1,19 | sort -u | cut -f2 | sort | uniq -c | sort -nrk1 > nr.interim.tsv

#Windows approach (Maybe not necessary)

python3 make_sliding_windows_all.py Ecoli.nr.filtered.prot.longestrepresentative.reduced2.faa --win 50 --step 10 > Ecoli.nr.filtered.prot.longestrepresentative.reduced2.windows.faa

#Python code:

#!/usr/bin/env python3
"""
make_sliding_windows.py
Generate sliding window fragments for *all* protein sequences in a FASTA file.

Usage:
  python3 make_sliding_windows.py all_proteins.faa --win 60 --step 10 > all_windows.faa
"""

import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Generate sliding windows for all sequences.")
parser.add_argument("fasta", help="Input FASTA file (protein sequences)")
parser.add_argument("--win", type=int, default=60, help="Window size (aa)")
parser.add_argument("--step", type=int, default=10, help="Step size (aa)")
args = parser.parse_args()

for record in SeqIO.parse(args.fasta, "fasta"):
    seq = str(record.seq)
    seq_len = len(seq)
    if seq_len < args.win:
        continue  # skip too-short sequences
    for i in range(0, seq_len - args.win + 1, args.step):
        frag = seq[i:i + args.win]
        # format: >originalID_start_end
        sys.stdout.write(f">{record.id}_{i+1}_{i+args.win}\n{frag}\n")

#Search windows against all proteins: 

/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Ecoli.nr.filtered.prot.longestrepresentative.reduced2.faa --db Ecoli.nr.filtered.prot.longestrepresentative.reduced2
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli.nr.filtered.prot.longestrepresentative.reduced2.windows.faa -d Ecoli.nr.filtered.prot.longestrepresentative.reduced2 --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore --ultra-sensitive --out Ecoli.nr.filtered.prot.longestrepresentative.reduced2.windows_vs_all.tsv -k 0 -b8 -c1
grep "^>" Ecoli.nr.filtered.prot.longestrepresentative.reduced2.faa | tr -d ">" > Ecoli.nr.filtered.prot.longestrepresentative.reduced2.txt
awk 'NR==FNR{keep[$1]=1; next} ($5>70 && $3>50) {n=split($1,a,"_"); id=a[1]; for(i=2;i<=n-2;i++) id=id"_"a[i]; if((id in keep) && id!=$2) print id "\t" $2}' Ecoli.nr.filtered.prot.longestrepresentative.reduced2.txt Ecoli.nr.filtered.prot.longestrepresentative.reduced2.windows_vs_all.tsv > windows_interim1.tsv

awk '{if ($1 < $2) print $1, $2; else print $2, $1}' windows_interim1.tsv | sort -u | sed "s/ /\t/g" > windows_interim2.tsv


