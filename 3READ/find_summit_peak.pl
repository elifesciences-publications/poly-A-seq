#！usr/bin/perl -w
use List::Util qw(reduce max);
use List::MoreUtils qw(uniq);
use strict;
use Getopt::Long;
use List::Compare;
#use Cwd 'abs_path';

my ($sam, $genome, $refflat, $utr_range, $bin, $threshold);

$sam = "";
$genome = "";
$refflat = "";
$utr_range = 1000;
$bin = 50;     # this number defines the left and right length that flank the cleavage site
$threshold = 1;

GetOptions 
(
  's=s'=>\$sam,		# name of sample
  'g=s'=>\$genome,	# path and file name of genome
  'ref=s'=>\$refflat,		# refflat file, can be generated from gtf or gff 
    );

if ($sam eq "" or $genome eq "" or $refflat eq "") {
	die "Missing one or more input file name(s)!
	use of script:
	perl findPAS.pl -s <sample_name> -g <path to genome fasta file ) -ref <refflat file>\n"} # check parameters

#$file_in =~ s/^.+\///g;
#$refflat =~ s/^.+\///g;

# read file from wig
my %plus; my %minus;
open (hand1, "$sam/thout/pos_plus.bedgraph") or die $!;
while (<hand1>)   {
	chomp;
	my ($chr, $start, $end, $val) = split /\t/;
	for my $i ($start..$end-1)    {
		$plus{$chr}{$i} = $val;   } }
close hand1;

open (hand1, "$sam/thout/pos_neg.bedgraph") or die $!;
while (<hand1>)   {
	chomp;
	my ($chr, $start, $end, $val) = split /\t/;
	for my $i ($start..$end-1)    {
		$minus{$chr}{$i} = $val; }  }
close hand1;

my $ge; my %genome;
open (hand2, $genome) or die $!;
while (<hand2>)   {
	$_=~ s/\s+$//;
    if (/^>/)       {
    $ge=$_;
    $ge=~ s/^>//;
    next;}
    $_ =~ tr/atcg/ATCG/;
    $genome{$ge} .= $_; }

close hand2;

# read the warning file, which contains the annotation that are not conventional ORF
my %exc;
open (W, "ref/translate_warning.txt") or print "no warning file found!\n";
while (<W>)   {
	$_=~ s/\s+$//;
	if (/^>/)       {
		$_ =~ s/^>//;
		$exc{$_} = ""; }
	 }
close W;

print "## finish loading references\n";



## find the summit for each peak 
my %gc; my %peak;

open (hand1, "$sam/merged_pos.bed") or die $!;

while (<hand1>)   {
	chomp;
	next if /^seqname/;
	
	# check bed file, should only have 5 columns
	my @a = split /\t/;
	my $col = @a;
	if ($col != 5)  {
		die "the column number of your bed file is not as requested,
		the bed file should be generated with bedtools merge function\n"}
		
	my ($chr, $start, $end, $strand, $frq) = split /\t/;
	#print "$frq\n";
	my %v; my($p,$max);
	
	for my $i ($start..$end)  {
		if ($strand eq "+")  {
			($p,$max) = searchPeak($chr, $start, $end, \%plus);}
		elsif ($strand eq "-")  {
			($p,$max) = searchPeak($chr, $start, $end, \%minus);}
		}
	
	# IMPORTANT: avoid using hash in hash unless it is very necessary, because it would cause circular reference and comsume huge amount of RAM.	
	$peak{"$chr\t$p\t$strand"} = "$max\t$frq";
	}
	
# use refflat reference file the annotate the peaks, either at 3' end or internal CDS	
my %cds; my %utr; my %cds_intron;
open (REF, $refflat) or die $!;
open (hand2, ">$sam/polyA_CDS.txt");
open (hand3, ">$sam/polyA_UTR3.txt");
open (in_seq, ">$sam/polyA_CDS_intron.txt");
open (hand4, ">$sam/gc_CDS.txt");
open (hand5, ">$sam/gc_UTR.txt");
open (in_gc, ">$sam/gc_CDS_intron.txt");
open (hand6, ">$sam/CDS_UTR_PolyA_counts.txt");
print hand6 "Gene\tCDS_counts\tUTR_counts\tintron_counts\tCDS_size\tFull_size\tCDS_intron_size\tUTR5_counts\tUTR5_len\tUTR5_GC\n";

my $item = 0;
my $prog = 0;

while (<REF>)   {
   chomp;
   my @a = split /\t/;
   my @ins = split /,/, $a[9];
   my @ine = split /,/, $a[10];
   my $cds_count = 0;
   my $utr3_count = 0;
   my $utr5_count = 0;
   my $intron_count = 0;
   my $full = 0;
   my $CDS = 0;
   my $cds_intron = 0;
   my $utr3_seq;
   my $utr5_seq;

   # define CDS region, first get all exon
   my @cds;
    if ($a[7] - $a[6] > 0 and !exists $exc{"$a[0]\_$a[1]"} )   {
	   
		for my $i ( 0..$a[8]-1 )  {
			for my $j ($ins[$i]..$ine[$i]-1)  {
				push @cds, $j; }  }
		$full = @cds;
		my @utrL; my @utrR;
		@utrL = $a[4]..$a[6]-1 if $a[4] - $a[6] != 0;
		@utrR = $a[7]..$a[5] if $a[5] - $a[7] != 0;

		my $lcl = List::Compare->new(\@cds, \@utrL);
		@cds = $lcl->get_unique;
		# needs to sort by value, not sort by ASCII code
		@cds = sort {$a <=> $b} @cds;
		#print join(",",@cds)."\n";
		my $lcr = List::Compare->new(\@cds, \@utrR);
		@cds = $lcr->get_unique;
		@cds = sort {$a <=> $b} @cds;
		#print join(",",@cds)."\n";
		$CDS = @cds;
		$cds_intron = $a[7] - $a[6] - $CDS;
		# searching the coding region, including intron, and found whether the 
		# PAS region is within the coding region (i.e, the -30 to -10 is completely
		# inside coding region.
		foreach my $i ($a[6]..$a[7])  {
			if (exists $peak{"$a[2]\t$i\t$a[3]"})   {
				my($max,$frq) = split /\t/, $peak{"$a[2]\t$i\t$a[3]"};
								
				# check whether the PAS region of peak is within CDS
				my $pas_st = $i-30;
				my $pas_end =$i-10;
				if ($a[3] eq '-')  {
					$pas_st = $i+10;
					$pas_end = $i+30;}
				
				# not good for this: if ( grep( /^$pas_st$/, @cds) and grep( /^$pas_end$/, @cds) )  {
				# use list::compare (intersection)
				my @pas = $pas_st..$pas_end;
				my $lc1 = List::Compare->new(\@cds, \@pas);
				my @intersect = $lc1->get_intersection;
				
				if (scalar @intersect == scalar @pas)  {
					# sum up the reads within CDS
					$cds_count += $frq;
					
					#set condition for CDS peaks
					if ($frq >= $threshold)    { 	 
						my $seq = getSeq($a[2], $i, $a[3], \%genome);
						#print "CDS_$a[0]\t$i\t$frq\n";
						for my $i (0..length($seq)-1)  {
						   my $base = substr($seq, $i, 1);
						   $cds{$i}{$base} ++;  }
						print hand2 ">$a[0]\_$a[1]\_CDS\n$seq\n";  }
						    } 
				else {	# record the sequence that within intron
					$intron_count += $frq;
					if ($frq >= $threshold)    { 	 
						my $seq = getSeq($a[2], $i, $a[3], \%genome);
						#print "CDS_$a[0]\t$i\t$frq\n";
						for my $i (0..length($seq)-1)  {
						   my $base = substr($seq, $i, 1);
						   $cds_intron{$i}{$base} ++;  }
						print in_seq ">$a[0]\_$a[1]\_intron\n$seq\n";  }   
					}   }
		
		}      }
				   
	
   ##---------------------------------------------------------------------------------------                            
   ## this part is for 3' UTR region, search only the highest peak and report the sequence 
   # define 3' UTR region, extend the 3' end 1000 bp so as to include the some signal that outside the defined region. This is typically the truth since defined 3' end are often not the real end
   my $st_utr; my $end_utr; 
	if ($a[3] eq '-')     {
		$st_utr = $a[4] - $utr_range;
		$st_utr = 0 if $st_utr < 0;
		$end_utr = $a[6];}
	else   {
		$st_utr = $a[7];
		$end_utr = $a[5] + $utr_range;}

   $utr3_seq = substr($genome{$a[2]}, $st_utr, $end_utr-$st_utr);
   if ($a[3] eq '-') {  
		$utr3_seq = reverse $utr3_seq;
		$utr3_seq =~ tr/ATGC/TACG/;  }

   # define 5' UTR region
   my $st_utr5; my $end_utr5; 
	if ($a[3] eq '-')     {
		$st_utr5 = $a[5] + $utr_range;
		$st_utr5 = length($genome{$a[2]}) if $end_utr5 > $genome{$a[2]};
		$end_utr5 = $a[7];}
	else   {
		$st_utr5 = $a[4] - $utr_range;
		$st_utr5 = 0 if $st_utr5 < 0;
		$end_utr5 = $a[6];}

   $utr5_seq = substr($genome{$a[2]}, $st_utr5, $end_utr5-$st_utr5);
   if ($a[3] eq '-') {  
		$utr5_seq = reverse $utr5_seq;
		$utr5_seq =~ tr/ATGC/TACG/;  }

   # count the number of reads in 5' UTR region
   	for my $i ($st_utr5..$end_utr5)   {
		
		if ( exists $peak{"$a[2]\t$i\t$a[3]"} )   {
			my($max,$frq) = split /\t/, $peak{"$a[2]\t$i\t$a[3]"};
			$utr5_count += $frq;
						} }
	
	my %peak_utr; 
	# first put all peaks in UTR region into a array, then sort the array from high to low, 
	# then can pick the desired one, highest or second highest
	for my $i ($st_utr..$end_utr)   {
		
		if ( exists $peak{"$a[2]\t$i\t$a[3]"} )   {
			my($max,$frq) = split /\t/, $peak{"$a[2]\t$i\t$a[3]"};
			$utr3_count += $frq;
			
			# remove some peaks that have very low counts, then put in a hash for value ranking
			next if $frq < 10; 
			$peak_utr{$i} = $max;
			
						} }
			
	next if scalar keys %peak_utr == 0;
	# get the position of highest value 
	## method1, use reduce from List::Utils
	my $k1 = reduce { $peak_utr{$a} > $peak_utr{$b} ? $a : $b } keys %peak_utr;
	my $v1 = max values %peak_utr;
	## key based on the second highest value (NOTE: the position of $b and $a decided descending or ascending order
	#my $k2 = (sort { $peak_utr{$b} <=> $peak_utr{$a} } keys %peak_utr)[1];
	#print "p$k1=$v1; p$k2=$v2\n";
	
	my $seq = getSeq($a[2], $k1, $a[3], \%genome);
	for my $i (0..length($seq)-1)  {
		my $base = substr($seq, $i, 1);
		$utr{$i}{$base} ++;  }
			
	print hand3 ">$a[0]\_$a[1]\_UTR3\n$seq\n";

	my $nG = ($utr5_seq =~ tr/G//);
	my $nC = ($utr5_seq =~ tr/C//);
	my $utr5_len = length($utr5_seq);
	my $utr5GC = ($nG+$nC)/$utr5_len;
	
	if ( $a[7] - $a[6] > 0 and !exists $exc{"$a[0]\_$a[1]"} )   {
		print hand6 "$a[0]\_$a[1]\t$cds_count\t$utr3_count\t$intron_count\t$CDS\t$full\t$cds_intron\t$utr5_count\t$utr5_len\t$utr5GC\n";}
	
	$item ++;
	if ($item == 500)  {
		$prog += $item;		
		$item = 0;
		print "$prog processed\n";
			}

	     } 


close REF;
	
	
		 
close hand1; close hand2; close hand3;

my @nuc = ("A", "C", "G", "T");

print hand4 "pos\tA\tC\tG\tT\n";

foreach my $i (sort keys %cds)  {
	my $p = $i-$bin;
	print hand4 "$p";
	foreach my $j (@nuc) {
		if (exists $cds{$i}{$j})  {
			print hand4 "\t$cds{$i}{$j}";}
		else {
			print hand4 "\t0";}  }
	print hand4 "\n";  }
	
close hand4;




print hand5 "pos\tA\tC\tG\tT\n";

foreach my $i (sort keys %utr)  {
	my $p = $i-$bin;
	print hand5 "$p";
	foreach my $j (@nuc) {
		if (exists $utr{$i}{$j})  {
			print hand5 "\t$utr{$i}{$j}";}
		else {
			print hand5 "\t0";}  }
	print hand5 "\n";  }
	
close hand5;

print in_gc "pos\tA\tC\tG\tT\n";

foreach my $i (sort keys %cds_intron)  {
	my $p = $i-$bin;
	print in_gc "$p";
	foreach my $j (@nuc) {
		if (exists $cds_intron{$i}{$j})  {
			print in_gc "\t$cds_intron{$i}{$j}";}
		else {
			print in_gc "\t0";}  }
	print in_gc "\n";  }
	
close in_gc;


##==============================================================
## this part for subs





sub searchPeak     {
	my($chr, $start, $end, $ref) = @_;
	my %v; my @v;
	for my $i ($start..$end)    {
		if (exists $$ref{$chr}{$i})   {
				my $val = $$ref{$chr}{$i};
				push @v, $val;
				$v{$val} = $i; }  }
	my $max = max @v;
	my $pos = $v{$max}; 
	#print "$max\t$pos\n";
	return ($pos,$max);     }
				
				
sub getSeq     {
	my ($chr, $p, $strand, $ref) = @_;
	my $seq = substr($$ref{$chr}, $p-$bin, $bin*2+1);
	if ($strand eq "-")   {
		$seq = reverse $seq;
		$seq =~ tr/ATCG/TAGC/; }
	return $seq;   }	
	
	
	
	
	
	
	
