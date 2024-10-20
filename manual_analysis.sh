#For any de novo gene candidate, concatenate all of its blastn hit sequences into one file, then run macse on them.
#Do this separately for tblastn, fwiw

sed "s/>/>outside_/g" ../flanks_genus/AELFKIGO_03201*blastn_seq.faa > AELFKIGO_03201_all_blastn_seq.faa
sed "s/>/>fegusonii_/g" ../flanks_species/AELFKIGO_03201*fergusonii_blastn_seq.faa | head -n-4 >> AELFKIGO_03201_all_blastn_seq.faa
sed "s/>/>ruysiae_/g" ../flanks_species/AELFKIGO_03201*ruysiae_blastn_seq.faa | head -n-4 >> AELFKIGO_03201_all_blastn_seq.faa
sed "s/>/>marmotae_/g" ../flanks_species/AELFKIGO_03201*marmotae_blastn_seq.faa | head -n-4 >> AELFKIGO_03201_all_blastn_seq.faa
sed "s/>/>whittamii_/g" ../flanks_species/AELFKIGO_03201*whittamii_blastn_seq.faa | head -n-4 >> AELFKIGO_03201_all_blastn_seq.faa
sed "s/>/>albertii_/g" ../flanks_species/AELFKIGO_03201*albertii_blastn_seq.faa | head -n-4 >> AELFKIGO_03201_all_blastn_seq.faa
head -n-4 AELFKIGO_03201_blastn_seq.faa >> AELFKIGO_03201_all_blastn_seq.faa

ls AELFKIGO_03201_all_blastn_seq.faa | awk '{OFS=""}{print "java -jar \/stor\/work\/Ochman\/hassan\/tools\/macse_v2.06.jar -prog alignSequences -fs 1000 -fs_term 1000 -stop 100 -seq ",$1," -out_NT ",$1,"_NT.aln -out_AA ",$1,"_AA.aln"}' | bash
