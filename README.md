# poly-A-seq

These codes were written to analyze PCPA (prematural cleavage and polyadenylation) at the 3' UTR, intron and CDS. 
For details of the analyses, please refer to our paper.
The use of arguments were written in the each scripts and can be changed for your need. To replicate the results in our publication, please use the settings described in our paper.  

## Before start
** require Perl(5.22) and R(3.4.0) or above**

**1, download and install bowtie2, tophat2, SAMTOOLS, BEDTools.**

**2, create the reference files (ORF sequences, CAI and CBI values).** 

Please go to iGenome (Illumnia) or other sources to download your desired genomic references, including genome fasta file, GTF and refFlat.
*use the perl script to create ORF reference sequences.*
open the script _sequence_bu_refFlat.pl_, make sure to give it the right path to access reference files. The run the code in termial
```
perl sequence_by_refFlat.pl 
```
Then put script _CBI_ref_define.R_ and _cal_CBIAvg.pl_ inside the *ref* folder. excute the R scripts in Rstudio and Perl in terminal as showed above
*Note*: install necessary module for R and Perl used in the codes. 

**3, download rest scripts and place them in the folder that contains your fastq file. **

## step1 to step4， processing reads and mapping.  
change the file name according inside the scripts (Shell or Perl), with a text editor. Then execute these code as below:
```
perl step1_trim3.pl
sh step2_align.sh
perl step3_PASS.pl
sh step4_position.sh
```
The above codes were written for **3READ** sequencing method. If **2P-seq** are used, please use the codes
followed by *\_2P-seq.pl 
## step5
change the file name inside the shell script and run, which should produce figures and files in subfolders
```
sh step5_PAS_search.pl
```
## step6
excute as following to analyze dicodon frequencies
```
sh step6_PAS_dicodon.sh
```
then execute _step6_PAS_dicodon.R_ in Rstudio software.
## step7
simply execute the shell script
```
sh step7_PAS_prediction.sh
```
## step8
execute the following script to calculate the PAS scores based on PSSM
```
perl step9_PSSM.pl
```

**NOTE: need to install required modules or packages for Perl and R**
Please direct all your questions to Yunkun Dang (email: yunkun_dang@126.com)

Physical address: Center of Life Science, College of Biological Sciences, Yunnan University, Kunming, Yunnan, China. 
