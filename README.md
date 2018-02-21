# poly(A)-seq pipeline

These codes were written to analyze PCPA (prematural cleavage and polyadenylation) at the 3' UTR, intron and CDS. 
For details of the analyses, please refer to our paper titled "Codon usage biases co-evolve with the transcription termination machinery to suppress premature cleavage and polyadenylation in coding regions".

The use of arguments were written in the each scripts and can be changed for your need. To replicate the results in our publication, please use the settings described in "Materials and Methods".  

## Before start
**require Perl(5.22) and R(3.4.0) or above**

**1, download and install bowtie2, tophat2, SAMTOOLS, BEDTools.**

make sure the installed softwares are in your PATH.

**2, create the reference files (ORF sequences, CAI and CBI values).** 

Please go to iGenome (Illumnia) or other sources to download your desired genomic references, including genome fasta file, GTF and refFlat.

*use the perl script to create ORF reference sequences.*

open the script _sequence_bu_refFlat.pl_, make sure to give it the right path to access reference files. The run the code in termial
```
perl sequence_by_refFlat.pl 
```
Then put script _CBI_ref_define.R_ and _cal_CBIAvg.pl_ inside the *ref* folder. excute the R scripts in Rstudio and Perl in terminal as showed above to create CAI and CBI reference file.

**Note**: install necessary module for R and Perl used in the codes

For R: **ggplot2 reshape2**

For Perl: **List::MoreUtil List::Compare Math::CDF** 

**3, put all your fastq files in a folder named "fq". Then download all scripts from "3READ" folder and place them outside the "fq" folder**

## step1 to step4， processing reads and mapping.  
change the file name according inside the scripts (Shell or Perl), with a text editor. Then execute these code as below:
```
perl step1_trim3.pl
sh step2_align.sh
perl step3_PASS.pl
sh step4_position.sh
```
The above codes were written for **3READS** sequencing method. If **2P-seq** are used, please use the codes
followed by *\_2P-seq.pl 
## step5
This step will create:
1. scatter plots of CBI/CAI to normalized ORF/3'UTR termination ratio
2. line plots showing the GC content flanking the pA sites within ORFs or 3' UTRs
HOW TO USE: change the file name inside the shell script and run, which should produce figures and files in subfolders
```
sh step5_PAS_search.pl
```
## step6
This step will create: codon and dicodon frequency around the putative PAS region.  
HOW TO USE: open the shell script and change the file names accordingly. Then excute as following to analyze codon and codon frequencies
```
sh step6_PAS_dicodon.sh
```
then execute _step6_PAS_dicodon.R_ in Rstudio software.
## step7
This step will create: 
1.The line plot of AU content flanking selected 6nt PAS motifs.
2. relative codon refrequency flanking the PAS motifs. 
HOW TO USE: Open and change the file name accordingly. Then simply execute the shell script, which will produce the figures 
```
sh step7_PAS_prediction.sh
```
## step8
This step will create the PSSM scores for PAS region of 3'UTR and ORF regions. 
execute the following script to calculate the PAS scores based on PSSM
```
perl step8_PSSM.pl
```
To make the boxplot, simly use the built-in function "boxplot" in R. 

***
Please direct all your questions to Yunkun Dang (email: yunkun_dang@126.com)

Physical address: Center for Life Science, College of Biological Sciences, Yunnan University, Kunming, Yunnan, China. 
