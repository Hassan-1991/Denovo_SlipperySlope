#Step-1:

for ani in $(seq 99.99 -0.1 96.0); do
    echo "Processing ANI cutoff $ani ..."

    awk -v cutoff="$ani" '$3 >= cutoff {print $1, $2}' E.coli_pairwise_ANI.tsv \
        > E.coli_pairwise_ANI.${ani}.pairs.txt

    python3 cluster_ani.py E.coli_pairwise_ANI.${ani}.pairs.txt \
        > E.coli_pairwise_ANI.${ani}.clusters.txt

    cat -n E.coli_pairwise_ANI.${ani}.clusters.txt \
        | sed "s/^ *//g" \
        | sed "s/ /\t/g" \
        | awk '{id=$1; for(i=2;i<=NF;i++) print id, $i}' \
        > E.coli_pairwise_ANI.${ani}.clusters

    total=$(wc -l < genome_names.tsv)
    cut -f2 -d " " E.coli_pairwise_ANI.${ani}.clusters \
        | grep -v -w -F -f - genome_names.tsv \
        | awk -v n="$total" '{print n--, $0}' - \
        >> E.coli_pairwise_ANI.${ani}.clusters
done

#Python code:

#!/usr/bin/env python3
import sys
from collections import defaultdict, deque

edges = defaultdict(set)
for line in open(sys.argv[1]):
    a, b = line.strip().split()
    edges[a].add(b)
    edges[b].add(a)

visited = set()
for node in edges:
    if node not in visited:
        cluster = []
        queue = deque([node])
        while queue:
            n = queue.popleft()
            if n not in visited:
                visited.add(n)
                cluster.append(n)
                queue.extend(edges[n])
        print(" ".join(sorted(cluster)))


######

#Count the number of genomes per cluster at different ANIs, deduplicated at 99.99%

sort -k2 E.coli_pairwise_ANI.99.99.clusters > interim
for i in $(ls *clusters | rev | cut -f2,3 -d "." | rev)
do
sort -k2 E.coli_pairwise_ANI."$i".clusters | join -1 2 -2 2 - interim | cut -f2- -d " " | sort -u | awk '{print $1}' | sort | uniq -c | awk '($1!=1)' | awk '{print "cluster"$2"\t"$1}' | sed "s/^/ANI"$i"\t/g" >> E.coli_pairwise_ANI.clusterinfo.dedup.tsv
sort -k2 E.coli_pairwise_ANI."$i".clusters | join -1 2 -2 2 - interim | cut -f2- -d " " | sort -u | awk '{print $1}' | sort | uniq -c | awk '($1==1)' | wc -l | sed "s/^/cluster_singleton\t/g" | sed "s/^/ANI"$i"\t/g" >> E.coli_pairwise_ANI.clusterinfo.dedup.tsv
done

#Now use this data to generate some images in R

#!/usr/bin/env Rscript

# Plot ANI cluster summaries
# Input: E.coli_pairwise_ANI.clusterinfo.dedup.tsv (cutoff, cluster, count)
# Output: side-by-side figure:
#   (left) stacked fractions per cutoff (singleton in red)
#   (right) mirrored bars: #clusters>=10 (top) and %genomes in clusters>=10 (bottom)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(tidyr)
  library(patchwork)
})

# ---- user params ----
infile <- "E.coli_pairwise_ANI.clusterinfo.dedup.tsv"
total_genomes <- 3331
# ---------------------

df <- read_tsv(
  infile,
  col_names = c("cutoff", "cluster", "count"),
  show_col_types = FALSE
) %>%
  mutate(
    cutoff = as.factor(cutoff),
    cluster = as.character(cluster),
    fraction = count / total_genomes
  ) %>%
  group_by(cutoff) %>%
  arrange(desc(fraction), .by_group = TRUE) %>%
  mutate(order_rank = row_number()) %>%
  ungroup() %>%
  mutate(order_flag = if_else(cluster == "cluster_singleton", Inf, order_rank)) %>%
  group_by(cutoff) %>%
  arrange(order_flag, .by_group = TRUE) %>%
  mutate(cluster = fct_reorder(cluster, fraction, .fun = identity)) %>%
  ungroup()

# colors: singleton red, everything else grey
clusters <- unique(as.character(df$cluster))
colors <- setNames(rep("grey70", length(clusters)), clusters)
if ("cluster_singleton" %in% names(colors)) colors["cluster_singleton"] <- "red"

p1 <- ggplot(df, aes(x = cutoff, y = fraction, fill = cluster)) +
  geom_col(width = 0.8, color = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = colors, guide = "none") +
  labs(x = "ANI cutoff", y = paste0("Fraction of genomes (n = ", total_genomes, ")")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )

df10 <- df %>%
  filter(cluster != "cluster_singleton") %>%
  group_by(cutoff) %>%
  summarise(
    n_clusters_ge10 = sum(count >= 10),
    frac_genomes_in_ge10 = sum(count[count >= 10]) / total_genomes,
    .groups = "drop"
  ) %>%
  mutate(frac_pct_neg = -(frac_genomes_in_ge10 * 100)) %>%
  pivot_longer(
    cols = c(n_clusters_ge10, frac_pct_neg),
    names_to = "metric",
    values_to = "value"
  )

p2 <- ggplot(df10, aes(x = cutoff, y = value, fill = metric)) +
  geom_col(width = 0.8) +
  scale_y_continuous(
    limits = c(-100, 55),
    breaks = c(seq(-100, 0, 20), seq(0, 55, 10)),
    labels = function(x) abs(x)
  ) +
  scale_fill_discrete(guide = "none") +
  labs(
    x = "ANI cutoff",
    y = "# clusters ≥10 (top) | % genomes in clusters ≥10 (bottom)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )

(p1 | p2)

#These figures tell me at which ANI:

#1 - Singletons are low
#Large (>=10 genome) clusters are high
#Large clusters cover a large number of total genomes

#ANI 99.49 has the highest number of large clusters (51), Lower range of singletons (8%), and high range of total pangenome coverage (~76.5%)
