python3 make_sliding_windows_all.py Ecoli.nr.filtered.prot.longestrepresentative.prot.faa --win 50 --step 10 > Ecoli.nr.filtered.prot.longestrepresentative.prot.windows.faa

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


