#!bin/sh

species='mouse_mm10'
len=50
start=30
down=50
dn=10

for sam in "WT" "Paf1" "Cdc73" "Ski8"

#"SRR2225336" "SRR2225337" "SRR2225342" "SRR2225343" "SRR2225338" "SRR2225339" "SRR2225340" "SRR2225341"  "SRR2225344" "SRR2225345" "SRR2225346" "SRR2225347" "SRR2225348" "SRR2225349" "SRR2225350" "SRR2225351"

	do
	echo "=================== processing $sam ======================="
	#cd $sam/thout
	#rm $sam\_filtered.bed $sam\_filtered_pos.bed $sam\_uniq.bed $sam\_filtered_pos_sorted.bed
	#cd ..
	#cd ..
	mkdir $sam/classify/

	echo "\t## merge reads into region"
	#bedtools merge -i $sam/thout/$sam\_sorted_pos.bed -s -d 24 -c 4 -o count > $sam/merged_pos.bed

	#perl 7.1_PAS_position.pl -s $sam -g /home/yunkun/genomes/$species/Sequence/WholeGenomeFasta/genome.fa -ref $sam/polyA_UTR3_-30--10_sig_motifs.txt -tr 2
	# tr is the number of motifs for analyses. 
	# this step is to generate the position of putative PAS motifs on genome


	#perl 7.2_PAS_classify.pl -s $sam -g /home/yunkun/genomes/$species/Sequence/WholeGenomeFasta/genome.fa -ref /home/yunkun/genomes/$species/Genes/refFlat.txt -tr 1 -len 80 -dn 80 

	Rscript 7.3_drawlines.R $sam 

	#perl 7.4_cleavage_classify.pl -s $sam -g /home/yunkun/genomes/$species/Sequence/WholeGenomeFasta/genome.fa -ref /home/yunkun/genomes/$species/Genes/refFlat.txt

        #mkdir $sam/pictures
	#cd $sam/pictures

	#kpLogo ../PAS_codons.fa -o PAS_codon
	#kpLogo ../PAS_amino_acid.fa -alphabet protein -o PAS_AA
	#kpLogo ../polyA_UTR3.txt -o PAS_UTR
	
	#cd ..
	#cd ..
	
	done
