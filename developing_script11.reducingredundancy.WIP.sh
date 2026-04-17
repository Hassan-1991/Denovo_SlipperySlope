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


