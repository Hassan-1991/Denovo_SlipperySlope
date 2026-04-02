cat all_genes_of_interest.noncoding.tsv ../all_genes_of_interest.presence.tsv  | cut -f4 | grep -v "Ecoli@" | sort -u > extraspecies.txt

cat all_genes_of_interest.noncoding.tsv ../all_genes_of_interest.presence.tsv | awk '($4!~"Ecoli@")' > all_gene_status.extraspecies.tsv
cat all_genes_of_interest.noncoding.tsv ../all_genes_of_interest.presence.tsv | awk '($4~"Ecoli@")' > all_gene_status.pangenome.tsv

#Compile:

awk -F '\t' 'BEGIN{OFS="\t"}
file==1 {
    genes[$1]
    next
}
file==2 {
    taxa[$1]
    next
}
file==3 && ($1 in genes) {
    status=$3
    if(status=="partial" || status=="ORFpartial") status="noncoding"
    else if(status=="ORFpresent") status="present"

    key=$1 SUBSEP $4

    if(!(key in best)) {
        best[key]=status
    } else if(status=="present") {
        best[key]="present"
    } else if(status=="noncoding" && best[key]!="present") {
        best[key]="noncoding"
    } else if(status=="flankpresent" && best[key]!="present" && best[key]!="noncoding") {
        best[key]="flankpresent"
    }
}
END{
    for(g in genes) {
        for(t in taxa) {
            key=g SUBSEP t
            if(key in best) print g, t, best[key]
            else print g, t, "absent"
        }
    }
}' file=1 ../step1_ORFans.txt file=2 extraspecies.txt file=3 all_gene_status.extraspecies.tsv > all_gene_status.extraspecies.compiled.tsv

awk -F '\t' 'BEGIN{OFS="\t"}
file==1 {
    genes[$1]
    next
}
file==2 {
    genomes[$1]
    next
}
file==3 {
    if($2 in genomes && !( $2 in taxon )) taxon[$2]=$4

    if($1 in genes) {
        status=$3
        if(status=="partial" || status=="ORFpartial") status="noncoding"
        else if(status=="ORFpresent") status="present"

        key=$1 SUBSEP $2

        if(!(key in best)) {
            best[key]=status
        } else if(status=="present") {
            best[key]="present"
        } else if(status=="noncoding" && best[key]!="present") {
            best[key]="noncoding"
        } else if(status=="flankpresent" && best[key]!="present" && best[key]!="noncoding") {
            best[key]="flankpresent"
        }
    }
}
END{
    for(g in genes) {
        for(genome in genomes) {
            key=g SUBSEP genome
            if(key in best) print g, genome, taxon[genome], best[key]
            else print g, genome, taxon[genome], "absent"
        }
    }
}' file=1 ../step1_ORFans.txt \
   file=2 <(cut -f1 ../../genome_cluster.contigs.tsv | sort -u) \
   file=3 all_gene_status.pangenome.tsv \
> all_gene_status.pangenome.compiled.tsv

#First transformation step:

gawk -F '\t' 'BEGIN{OFS=""}
FNR==NR {
    lineage[++nlineage] = $1
    next
}
{
    gene = $1
    genome = $2
    taxon = $3
    stat = $4

    if (!(gene in seen_gene)) {
        seen_gene[gene] = 1
        genes[++ngenes] = gene
    }

    status[gene SUBSEP genome] = stat
    taxgen[taxon SUBSEP genome] = 1
}
END {
    for (i=1; i<=nlineage; i++) {
        tax = lineage[i]
        delete tmp
        for (k in taxgen) {
            split(k, a, SUBSEP)
            if (a[1] == tax) tmp[a[2]] = 1
        }
        ngen[tax] = asorti(tmp, ord)
        for (j=1; j<=ngen[tax]; j++) {
            genome_order[tax, j] = ord[j]
        }
    }

    for (g=1; g<=ngenes; g++) {
        gene = genes[g]
        printf "%s", gene

        for (i=1; i<=nlineage; i++) {
            tax = lineage[i]
            printf ",%s(", tax

            for (j=1; j<=ngen[tax]; j++) {
                genome = genome_order[tax, j]
                key = gene SUBSEP genome
                printf "%s%s", (j>1 ? "," : ""), (key in status ? status[key] : "absent")
            }

            printf ")"
        }

        printf "\n"
    }
}' lineage_order.txt all_gene_status.pangenome.compiled.tsv > all_gene_status.pangenome.nested.csv

#Second transformation step (Python):

#!/usr/bin/env python3

import re
from collections import OrderedDict

INPUT = "all_gene_status.pangenome.nested.csv"
RED_OUT = "allredcounts.csv"
GREEN_OUT = "allgreencounts.csv"
YELLOW_OUT = "allyellowcounts.csv"
MAP_OUT = "gene_to_distribution_map.tsv"


def split_top_level_commas(s):
    parts = []
    buf = []
    depth = 0
    for ch in s.rstrip("\n"):
        if ch == "," and depth == 0:
            parts.append("".join(buf))
            buf = []
            continue
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        buf.append(ch)
    parts.append("".join(buf))
    return parts


def parse_block(block):
    m = re.match(r"^(.*?)\((.*)\)$", block)
    if not m:
        raise ValueError(f"Could not parse block: {block}")
    taxon = m.group(1)
    statuses = m.group(2).split(",") if m.group(2) else []
    return taxon, statuses


records = []

with open(INPUT) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue

        fields = split_top_level_commas(line)
        gene = fields[0]
        blocks = fields[1:]

        red = []
        green = []
        yellow = []

        total_present = 0

        for block in blocks:
            taxon, statuses = parse_block(block)

            r = sum(1 for x in statuses if x == "noncoding")
            g = sum(1 for x in statuses if x == "present")
            y = sum(1 for x in statuses if x == "flankpresent")

            red.append(str(r))
            green.append(str(g))
            yellow.append(str(y))

            total_present += g

        red_row = ",".join(red)
        green_row = ",".join(green)
        yellow_row = ",".join(yellow)

        dist_key = red_row + "||" + green_row + "||" + yellow_row

        records.append({
            "gene": gene,
            "total_present": total_present,
            "red": red_row,
            "green": green_row,
            "yellow": yellow_row,
            "dist_key": dist_key
        })

# sort by total present descending, then gene name for stability
records.sort(key=lambda x: (-x["total_present"], x["gene"]))

# assign unique distribution IDs preserving first appearance after sorting
unique_rows = OrderedDict()
for rec in records:
    if rec["dist_key"] not in unique_rows:
        unique_rows[rec["dist_key"]] = {
            "distribution_id": f"D{len(unique_rows)+1}",
            "red": rec["red"],
            "green": rec["green"],
            "yellow": rec["yellow"],
            "representative_gene": rec["gene"]
        }

# write mapping file
with open(MAP_OUT, "w") as out:
    out.write("gene\ttotal_present\tdistribution_id\tred_row\tgreen_row\tyellow_row\n")
    for rec in records:
        dist_id = unique_rows[rec["dist_key"]]["distribution_id"]
        out.write(
            f"{rec['gene']}\t{rec['total_present']}\t{dist_id}\t"
            f"{rec['red']}\t{rec['green']}\t{rec['yellow']}\n"
        )

# write deduplicated count matrices
with open(RED_OUT, "w") as redf, open(GREEN_OUT, "w") as greenf, open(YELLOW_OUT, "w") as yellowf:
    for item in unique_rows.values():
        redf.write(item["red"] + "\n")
        greenf.write(item["green"] + "\n")
        yellowf.write(item["yellow"] + "\n")

python3 make_count_matrices.py

#Finally, the figure plot. This is run on google colab, but could be run on terminal as well:

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

max_dots = 10

# Load counts
red_counts = pd.read_csv('/content/allredcounts.csv', header=None).values.T
green_counts = pd.read_csv('/content/allgreencounts.csv', header=None).values.T
yellow_counts = pd.read_csv('/content/allyellowcounts.csv', header=None).values.T

# Infer grid size from file
rows, cols = red_counts.shape

# Sanity checks
assert green_counts.shape == (rows, cols)
assert yellow_counts.shape == (rows, cols)

# Make sure no cell exceeds 10 total dots
total_counts = red_counts + green_counts + yellow_counts
mask = total_counts > max_dots

# Adjust yellow downward first, then green if needed
yellow_counts[mask] = np.maximum(0, max_dots - red_counts[mask] - green_counts[mask])
green_counts[mask]  = np.maximum(0, max_dots - red_counts[mask] - yellow_counts[mask])

# Recompute total to confirm
total_counts = red_counts + green_counts + yellow_counts
assert np.all(total_counts <= max_dots)

# Plot
fig, ax = plt.subplots(figsize=(cols, rows))

box_size = 1
ball_radius = 0.15

offsets = [
    (-0.25, 0.25), (-0.25, 0), (-0.25, -0.25),
    (0, 0.25),     (0, 0),     (0, -0.25),
    (0.25, 0.25),  (0.25, 0),  (0.25, -0.25)
]

for i in range(rows):
    for j in range(cols):

        red = red_counts[i, j]
        green = green_counts[i, j]
        yellow = yellow_counts[i, j]
        total = red + green + yellow

        if total == max_dots:
            if red >= green and red >= yellow:
                facecolor = 'lightcoral'
            elif green >= red and green >= yellow:
                facecolor = 'lightgreen'
            else:
                facecolor = 'khaki'
        else:
            facecolor = 'white'

        rect = patches.Rectangle(
            (j, rows - i - 1),
            box_size, box_size,
            linewidth=0.5,
            edgecolor='black',
            facecolor=facecolor
        )
        ax.add_patch(rect)

        if total < max_dots:
            for k in range(red):
                dx, dy = offsets[k] if k < 9 else (0, 0)
                ax.add_patch(
                    patches.Circle((j + 0.5 + dx, rows - i - 0.5 + dy),
                                   radius=ball_radius, color='red')
                )

            for k in range(green):
                idx = red + k
                dx, dy = offsets[idx] if idx < 9 else (0, 0)
                ax.add_patch(
                    patches.Circle((j + 0.5 + dx, rows - i - 0.5 + dy),
                                   radius=ball_radius, color='green')
                )

            for k in range(yellow):
                idx = red + green + k
                dx, dy = offsets[idx] if idx < 9 else (0, 0)
                ax.add_patch(
                    patches.Circle((j + 0.5 + dx, rows - i - 0.5 + dy),
                                   radius=ball_radius, color='gold')
                )

ax.set_xlim(0, cols)
ax.set_ylim(0, rows)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.show()

#Exproting from colab:

fig.set_size_inches(cols, rows)  # Set the size again to ensure big dimensions
fig.savefig("dot_matrix_plot.updated.111325.pdf", bbox_inches='tight')

from google.colab import files
files.download("dot_matrix_plot.updated.111325.pdf")
