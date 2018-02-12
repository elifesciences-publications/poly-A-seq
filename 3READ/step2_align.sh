#!/bin/sh


### created and update on 04-15-2014 by Yunkun as a ribosomal profiling pipeline


# Print start status message.
echo "job started"
start_time=`date +%s`

species="mouse_mm10"

#defile your file name (in shell, define variables do not use $)



for input in "SRR2225336" "SRR2225337" "SRR2225338" "SRR2225339" "SRR2225340" "SRR2225341" "SRR2225342" "SRR2225343" "SRR2225344" "SRR2225345" "SRR2225346" "SRR2225347" "SRR2225348" "SRR2225349" "SRR2225350" "SRR2225351"

do 

mkdir $input

tophat2 -p 8 -G ~/genomes/$species/Genes/genes.gtf -o $input/thout --library-type=fr-secondstrand ~/genomes/$species/Sequence/Bowtie2Index/genome fq/$input\_A_1.fq

# secondstrand use for ligation. The fq data has been reverse complemented
# strand information:  fr-firststrand	dUTP, NSR, NNSR	Same as above except we enforce the rule that the right-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during first strand synthesis is sequenced.
## here although the method itself is ligation based, but the sequencing result is reversed complement. So we use fr-firststrand


#echo "$input alignment finished"
cd $input
cd thout
samtools index accepted_hits.bam

# remove the reads that has multiple targets
samtools view -bq 1 accepted_hits.bam > $input\_uniq.bam

#samtools view -BS $input/thout/accepted_hits.bam | grep -E "@|NM:" | grep -v "XS:" > $input/thout/unique_mapper.sam

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
