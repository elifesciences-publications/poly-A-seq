#!bin/sh
# set the species used for analysis, in shell, no space between varible and = 
species='NC10'

for sam in "Nuc_1" "Nuc_2" "Total_1" "Total_2"

do 

# need to create a chrom file 

echo "===========processing $sam===========\n"
bedtools genomecov -i $sam/thout/$sam\_filtered_pos.bed -bg -g /home/yunkun/genomes/$species/Sequences/WholeGenomeFasta/genome.chrom -strand + > $sam/thout/pos_plus_p.bedgraph

bedtools genomecov -i $sam/thout/$sam\_filtered_pos.bed -bg -g /home/yunkun/genomes/$species/Sequences/WholeGenomeFasta/genome.chrom -strand - > $sam/thout/pos_neg_p.bedgraph

# use bedgraph data as reference to justify which one is good or bad. Critiria: 
# 1, look at the 20nt sequence downstream of cleavage site. if AAAAAA or pattern like A{3,}.A{2,}.A{2,}, then remove it
# 2, if within a 12nt window within that 20nt region, more than 8A were found, then remove

perl PAS_filter.pl -s $sam -strand n -species $species

perl PAS_filter.pl -s $sam -strand p -species $species


# combine filtered bed files
cat $sam/thout/plus_filtered_PAS.bed $sam/thout/neg_filtered_PAS.bed > $sam/thout/PAS_filtered_pos.bed
rm $sam/thout/plus_filtered_PAS.bed $sam/thout/neg_filtered_PAS.bed

# create combined bedgraph for igv 
awk '{$4="-"$4; print $0}' $sam/thout/pos_neg_filtered.bedgraph > $sam/thout/pos_neg_1.bedgraph
cat $sam/thout/pos_neg_1.bedgraph $sam/thout/pos_plus_filtered.bedgraph > $sam/thout/pos_all.bedgraph

sort -k1,1 -k2,2n $sam/thout/pos_all.bedgraph > $sam/$sam\_all.wig

rm $sam/thout/pos_all.bedgraph $sam/thout/pos_neg_1.bedgraph $sam/thout/pos_plus_p.bedgraph $sam/thout/pos_neg_p.bedgraph
done

