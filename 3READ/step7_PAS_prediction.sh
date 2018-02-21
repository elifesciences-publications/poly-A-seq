#!bin/sh

species='mouse_mm10'


for sam in "SRR2225336" "SRR2225337" "SRR2225342" "SRR2225343" 

	do
	echo "=================== processing $sam ======================="
	#cd $sam/thout
	#rm $sam\_filtered.bed $sam\_filtered_pos.bed $sam\_uniq.bed $sam\_filtered_pos_sorted.bed
	#cd ..
	#cd ..
	mkdir $sam/classify/

	echo "\t## merge reads into region"

	perl 7.1_PAS_position.pl -s $sam -g ~/genomes/$species/Sequence/WholeGenomeFasta/genome.fa -ref $sam/polyA_UTR3_-30--10_sig_motifs.txt -tr 2
	# tr is the number of motifs for analyses. 
	# this step is to generate the position of putative PAS motifs on genome


	perl 7.2_PAS_classify.pl -s $sam -g /home/yunkun/genomes/$species/Sequence/WholeGenomeFasta/genome.fa -ref /home/yunkun/genomes/$species/Genes/refFlat.txt -tr 1 -len 80 -dn 80 

	Rscript 7.3_drawlines.R $sam 

   
	
	done
