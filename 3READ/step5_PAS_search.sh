#!bin/sh

# DEFINE the region where PAS motifs are supposed to locate. typically -30 to -10 upstream of cleavage site (pA site)
st="-30"
ed="-10"


species='mouse_mm10'

for sam in  "SRR2225336" "SRR2225337" "SRR2225338" "SRR2225339" "SRR2225340" "SRR2225341" "SRR2225342" "SRR2225343" "SRR2225344" "SRR2225345" "SRR2225346" "SRR2225347" "SRR2225348" "SRR2225349" "SRR2225350" "SRR2225351"

	do 
	echo "==================================\nnow processing $sam\n=================================="
	cd $sam
	# sort the bed file and merge into clusters with ones that within 24nt distance.
	# also count the frequency of peaks
	echo "\t## sorting bed files"
	sort -k1,1 -k2,2n thout/$sam\_filtered_pos.bed > thout/$sam\_filtered_pos_sorted.bed
	
	## -c the column number that is used for operation (-o), here the column number is started
	## from 1, not 0
	echo "\t## merge reads into region"
	bedtools merge -i thout/$sam\_filtered_pos_sorted.bed -s -d 24 -c 4 -o count > merged_pos.bed

	cd ..

	## the following scripts are for searching peaks and motifs
	echo "\t## searching peaks at UTR and CDS with perl script" 
	perl find_summit_peak.pl -s $sam -g <PATH_to_genome_fasta_file> -ref <PATH_to_refFlat_file>
	
	echo "\t## searching PAS signal in UTR"
	perl findPAS.pl -in $sam/polyA_UTR3.txt -start $st -end $ed
	
	echo "\t## searching PAS signal in CDS"
	perl findPAS.pl -in $sam/polyA_CDS.txt -start $st -end $ed

	echo "\t## searching PAS signal in intron in CDS region"
	perl findPAS.pl -in $sam/polyA_CDS_intron.txt -start $st -end $ed

	echo "\t## running R scripts to draw graphs"
	Rscript plot_polyA.R $sam <name_of_CAI_CBI_reference> $st $ed
	
	done 

