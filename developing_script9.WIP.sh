#First, let's compile both sets of sequences together

mkdir ORFpresence
cd ORFpresence
#Make a list of query proteins for which we have either flank-guided or unguided DNA sequence hits

for i in $(cat ../step1_ORFans.txt)
do
cat ../flank/"$i"_interval.fna | sed 's/:.*//' | sed "s/>/>flank_/g" >> "$i"_compiled_redundant.fna
cat ../regular/"$i"_blastn.seq.faa | sed "s/>/>regular_/g" >> "$i"_compiled_redundant.fna
done

#Now to extract ORFs from these sequences and search the query against them
for i in $(ls *_compiled_redundant.fna | rev | cut -f3- -d "_" | rev)
do
    echo "seqkit fx2tab \"${i}_compiled_redundant.fna\" | sed \"s/\t\$//g\" | sed \"s/^/>/g\" | sed \"s/\t/\n/g\" | awk '/^>/ {print} /^[^>]/ {gsub(/-/, \"\"); print}' | head -n-2 > \"${i}.ORFsearch_targets.frames.fna\"" >> "${i}.sh"
    echo "seqkit fx2tab \"${i}_compiled_redundant.fna\" | sed \"s/\t\$//g\" | sed \"s/^/>/g\" | sed \"s/\t/\n/g\" | awk '/^>/ {print} /^[^>]/ {gsub(/-/, \"\"); print}' | head -n-2 | awk '/^>/ {print} /^[^>]/ {print substr(\$0, 2)}' | sed \"s/>/>frame1_/g\" >> \"${i}.ORFsearch_targets.frames.fna\"" >> "${i}.sh"
    echo "seqkit fx2tab \"${i}_compiled_redundant.fna\" | sed \"s/\t\$//g\" | sed \"s/^/>/g\" | sed \"s/\t/\n/g\" | awk '/^>/ {print} /^[^>]/ {gsub(/-/, \"\"); print}' | head -n-2 | awk '/^>/ {print} /^[^>]/ {print substr(\$0, 3)}' | sed \"s/>/>frame2_/g\" >> \"${i}.ORFsearch_targets.frames.fna\"" >> "${i}.sh"
    echo "awk '/^>/ {print} /^[^>]/ {gsub(/.{3}/, \"& \"); print}' \"${i}.ORFsearch_targets.frames.fna\" | sed \"s/ATG/M/g\" | sed \"s/TTG/L/g\" | sed \"s/GTG/V/g\" | sed \"s/TAA/Z/g\" | sed \"s/TGA/Z/g\" | sed \"s/TAG/Z/g\" | sed \"s/ //g\" | awk '/^>/ {print} /^[^>]/ {print \$0 \"Z\"}' | awk '/^>/ {header = \$0} /^[^>]/ { while(match(\$0, /[MLV][^Z]*Z/)) { orf = substr(\$0, RSTART, RLENGTH); print header \": \" orf; \$0 = substr(\$0, RSTART + RLENGTH) } }' | cat -n | sed \"s/^ *//g\" | sed \"s/\t>/_ORF_/g\" | sed \"s/^/>/g\" | sed \"s/: /\n/g\" | awk '/^>/ {print} /^[^>]/ {gsub(\"M\", \"ATG\"); gsub(\"L\", \"TTG\"); gsub(\"V\", \"GTG\"); print}' > \"${i}.ORFsearch_targets.ORFs.fna\"" >> "${i}.sh"
    echo "/stor/work/Ochman/hassan/tools/faTrans \"${i}.ORFsearch_targets.ORFs.fna\" \"${i}.ORFsearch_targets.ORFs.prot.faa\"" >> "${i}.sh"
    echo "grep --no-group-separator -A1 -w \"${i}\" ../step1_ORFans.prot.faa > \"${i}.query.prot.faa\"" >> "${i}.sh"
    echo "/stor/work/Ochman/hassan/tools/fasta-36.3.8i/bin/ssearch36 -m 8 \"${i}.query.prot.faa\" \"${i}.ORFsearch_targets.ORFs.prot.faa\" > \"${i}.results.ssearch.tab\"" >> "${i}.sh"
done #parallelize however

#Process the results to get a list of all genomes (or taxa for extra-pangenome) which has a matching ORF, either full or partial
for i in $(ls *_compiled_redundant.fna | rev | cut -f3- -d "_" | rev)
do
length=$(seqkit fx2tab "$i".query.prot.faa | awk '{print length($2)}')
awk -F '\t' -v var="$length" '{cov = ($8 - $7 + 1) / var; print $0 "\t" cov}' "$i".results.ssearch.tab | awk -F '\t' '($NF>=0.6&&$11<0.001)' | cut -f2 | sed "s/frame1_//g" | sed "s/frame2_//g" | cut -f3- -d "_" | sed "s/_/\t/1" | sed "s/^/"$i"\t/g" >> ORF_presence.eval001.tsv
awk -F '\t' -v var="$length" '{cov = ($8 - $7 + 1) / var; print $0 "\t" cov}' "$i".results.ssearch.tab | awk -F '\t' '($NF<0.6&&$11<0.001)' | cut -f2 | sed "s/frame1_//g" | sed "s/frame2_//g" | cut -f3- -d "_" | sed "s/_/\t/1" | sed "s/^/"$i"\t/g" >> ORF_partial.eval001.tsv
done

#These files will be required eventually for purposes of compiling results

