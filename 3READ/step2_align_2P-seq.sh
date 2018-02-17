#!/bin/sh


### created and update on 04-15-2014 by Yunkun as a ribosomal profiling pipeline


# Print start status message.
echo "job started"
start_time=`date +%s`

species="NC_12"

#defile your file name (in shell, define variables do not use $)



for input in "CHX_1" "CHX_2" "Nuc_1" "Nuc_2" "Total_1" "Total_2"

do 

mkdir $input

tophat2 -p 8 -G ~/Deepseq/NC10_transcripts.gff3 -o $input/thout --library-type=fr-secondstrand ~/Deepseq/Indexes/NC10_b2 fq/$input\_trimA.fq

# secondstrand use for ligation. The fq data has been reverse complemented
# strand information:  fr-firststrand	dUTP, NSR, NNSR	Same as above except we enforce the rule that the right-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during first strand synthesis is sequenced.
## here although the method itself is ligation based, but the sequencing result is reversed complement. So we use fr-firststrand


#echo "$input alignment finished"
cd $input/thout

# remove the reads that has multiple targets
samtools view -bq 1 accepted_hits.bam > $input\_uniq.bam
samtools index $input\_uniq.bam

# create bed files
bedtools bamtobed -bed12 -i $input\_uniq.bam > $input\_uniq.bed

cd ..

cd ..

	done




# Print end status message.
echo
echo "job finished"
end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.
