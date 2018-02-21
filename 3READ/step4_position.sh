#!bin/sh

species='mouse_mm10'


# put your sample name here as used in step1. quote each name and seperate with space
for sam in "WT" "Paf1" "Cdc73" "Ski8"



do 

  echo "===========processing $sam===========\n"

  sort -k1,1 -k2,2n $sam/thout/$sam\_filtered_pos.bed > $sam/thout/$sam\_sorted_pos.bed
  
  # create chrom file if not exists. In short this file indicated the length of each chromosome. Please check BEDTools suite (genomecov) for more information.
  bedtools genomecov -i $sam/thout/$sam\_sorted_pos.bed -bg -g <chrom.file> -strand + > $sam/thout/pos_plus.bedgraph

  bedtools genomecov -i $sam/thout/$sam\_sorted_pos.bed -bg -g <chrom.file> -strand - > $sam/thout/pos_neg.bedgraph

  # create combined bedgraph for igv 
  awk '{$4="-"$4; print $0}' $sam/thout/pos_neg.bedgraph > $sam/thout/pos_neg_1.bedgraph
  cat $sam/thout/pos_neg_1.bedgraph $sam/thout/pos_plus.bedgraph > $sam/thout/pos_all.bedgraph

  sort -k1,1 -k2,2n $sam/thout/pos_all.bedgraph > $sam/$sam\_all.wig

  rm $sam/thout/pos_all.bedgraph $sam/thout/pos_neg_1.bedgraph $sam/thout/$sam\_filtered_pos.bed
  
done

