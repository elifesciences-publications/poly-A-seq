#!bin/sh

species='mouse_mm10'

for sam in "SRR2225336" "SRR2225337" "SRR2225338" "SRR2225339" "SRR2225340" "SRR2225341" "SRR2225342" "SRR2225343" "SRR2225344" "SRR2225345" "SRR2225346" "SRR2225347" "SRR2225348" "SRR2225349" "SRR2225350" "SRR2225351"


	do
	#cd $sam/thout
	#rm $sam\_filtered.bed $sam\_filtered_pos.bed $sam\_uniq.bed $sam\_filtered_pos_sorted.bed
	#cd ..
	#cd ..
	perl step6_PAS_codon.pl -s $sam -g /home/yunkun/genomes/$species/Sequence/WholeGenomeFasta/genome.fa -ref /home/yunkun/genomes/$species/Genes/refFlat.txt -start 45 -end -45

	
	done
	
	perl 6_count_dicodon.pl
