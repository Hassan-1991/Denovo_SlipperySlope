#All 500 this time
#Do clustering solely based on mmseqs2 criteria

cat ../500_gffs/500_proteins/*proteins.faa > all_500_proteins.faa
mmseqs createdb all_500_proteins.faa all_500_proteins
mmseqs search all_500_proteins all_500_proteins resultDB tmp --min-seq-id 0.8 -c 0.7 --cov-mode 1
mmseqs convertalis all_500_proteins all_500_proteins resultDB resultDB.m8
mmseqs linclust all_500_proteins clusterDB tmp --min-seq-id 0.8 -c 0.7 --cov-mode 1
mmseqs createtsv all_500_proteins all_500_proteins clusterDB all_500_proteins.clusters.tsv
