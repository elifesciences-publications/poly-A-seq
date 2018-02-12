#!usr/bin/perl -w
use List::Util qw(reduce max sum);
use List::MoreUtils qw(uniq);
use strict;
use Getopt::Long;
use List::Compare;
use Data::Dumper;
#use Cwd 'abs_path';
use Number::Format "round";



=head1
Purpose: 1, classify the prediction of PAS signal based on real 2P/3P seq data, i.e. which is the sequence causing transcription real stop or read though. 

=cut



my ($sam, $genome, $refflat, $reads, $length, $dn, $range_up, $range_down);

$sam = "";
$genome = "";
$refflat = "";
$reads = 1; # default number to select region (threshold for RNA)
$length = 50; # default length of upstream sequence to pick.
$dn=50; # default length of downstream sequence to pick, 0 is the point of polyA signal site.
$range_up = 5;	#here is the search range, from 10-40 after the the first nt of pattern. See whether cleavage signal is within this range. If yes, that means the PAS is likely a positive signal
$range_down = 35;

my $aa_up = 0;		# the codon number that upstream of 6nt PAS motif
my $aa_dn = 15;	# the codon number that downstream of 6nt PAS motif



GetOptions 
(
  's=s'=>\$sam,		# name of sample
  'g=s'=>\$genome,	# path and file name of genome
  'ref=s'=>\$refflat,	# refflat file, can be generated from gtf or gff 
  'tr=s'=>\$reads,
  'len=s'=>\$length,
  'dn=s'=>\$dn
  );

if ($sam eq "" or $genome eq "" or $refflat eq "") {
	die "Missing one or more input file name(s)!
	use of script:
	perl findPAS.pl -s <sample_name> -g <path to genome fasta file ) -ref <refflat file> \n"} # check parameters

#$file_in =~ s/^.+\///g;
#$refflat =~ s/^.+\///g;

print "++++ The parameter you choose ++++\n threshold=$reads, upstream seq length = $length \n";

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

my %orf;
open (hand3, "ref/CDS_DNA.fa") or die "no reference CDS DNA sequences file";
while (<hand3>)   {
	$_=~ s/\s+$//;
    	if (/^>/)       {
    		$ge=$_;
    		$ge=~ s/^>//;
    		next;}
	$_ =~ uc $_;
    	$orf{$ge} .= $_; }

close hand3;


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


## justify whether the downstread of a PAS may have cleavage signal or not

my (%true_pas, %false_pas);
open (hand1, "$sam/classify/polyA_prediction_plus.bed") or die $!;
while (<hand1>)   {
	chomp;
	my @a = split /\t/;
	$a[3] =~ s/_[+-]$//;
	my ($p,$max,$sum) = searchPeak($a[0], $a[2]+$range_up-1, $a[2]+$range_down-1, \%plus);
	if ($sum >= $reads)   {
		$true_pas{"$a[0]\t$a[1]\t+"} = "$a[3]\t$p\t$max\t$sum";
		}
	else {
		$false_pas{"$a[0]\t$a[1]\t+"} = $a[3];
		}
	}

close hand1;
	
open (hand1, "$sam/classify/polyA_prediction_minus.bed") or die $!;
while (<hand1>)   {
	chomp;
	my @a = split /\t/;
	$a[3] =~ s/_[+-]$//;
	my ($p,$max,$sum) = searchPeak($a[0], $a[1]-$range_down, $a[1]-$range_up, \%plus);
	#print "$a[0]\t$p\t$max\t$sum\n";
	if ($sum >= $reads)   {
		$true_pas{"$a[0]\t$a[1]\t-"} = "$a[3]\t$p\t$max\t$sum";
		}
	else {
		$false_pas{"$a[0]\t$a[1]\t-"} = $a[3];
		}
	}

close hand1;



# use refflat reference file the annotate the peaks, either at 3' end or internal CDS	
my %cdsT; my %cdsF; my %cds_intronT; my %cds_intronF; my %utrT; my %utrF;


open (REF, $refflat) or die $!;
open (trueCDS, ">$sam/classify/CDS_true.fa");
open (falseCDS, ">$sam/classify/CDS_false.fa");
open (trueIntron, ">$sam/classify/CDS_intron_true.fa");
open (falseIntron, ">$sam/classify/CDS_intron_false.fa");
open (trueUTR, ">$sam/classify/UTR_true.fa");
open (falseUTR, ">$sam/classify/UTR_false.fa");
open (codonT, ">$sam/classify/codon_$aa_up.UP_$aa_dn.Down_PAS_T.txt");
open (codonF, ">$sam/classify/codon_$aa_up.UP_$aa_dn.Down_PAS_F.txt");

my $item = 0;
my $prog = 0;
my $cds_countT = 0;
my $cds_countF = 0;
my $intron_countT = 0;
my $intron_countF = 0;
my $utr_countT = 0;
my $utr_countF = 0;



while (<REF>)   {
   chomp;
   my @a = split /\t/;
   my @ins = split /,/, $a[9];
   my @ine = split /,/, $a[10];


   # define coding region, first get all exon
   my %cr;  my %utr3;

	if ($a[7] - $a[6] > 0 and !exists $exc{"$a[0]\_$a[1]"} )   {  
	# coding region exists, i.e. not noncoding RNA
	   
		for my $i ( 0..$a[8]-1 )  {
			for my $j ($ins[$i]..$ine[$i]-1)  {
				$cr{$j} = "";
		 }  }
=head1		
		# remove utrL from hash cr
		if ($a[4] - $a[6] != 0)   {
			for my $i ($a[4]..$a[6]-1)  {
				delete $cr{$i} if exists $cr{$i};
		 } }

		
		# remove utrR from hash cr		
		if ($a[5] - $a[7] != 0)   {
			for my $i ($a[5]..$a[7]-1)  {
				delete $cr{$i} if exists $cr{$i};
		 } }
=cut
		
		}


   # work on coding region, including the ORF and intron within 

    if ($a[7] - $a[6] > 0 and !exists $exc{"$a[0]\_$a[1]"} )   {

		foreach my $i ($a[6]..$a[7])  {   # from ATG to stop codon, including introns

			####### true signal count #########
			if (exists $true_pas{"$a[2]\t$i\t$a[3]"})   {

				my($motif,$pos,$max,$frq) = split /\t/, $true_pas{"$a[2]\t$i\t$a[3]"};
				my $seq = getSeq($a[2], $i, $a[3], \%genome);							
				# check whether the PAS region of peak is within CDS
				my $pas_st = $i;
				my $pas_end = $i + 6;
				if ($a[3] eq '-')  {
					$pas_st = $i-6;
					$pas_end = $i;}

				my $overlap = 0;
				for my $j ($pas_st .. $pas_end)  {
					$overlap ++ if !exists $cr{$j}; }
				
				if ( $overlap == 0 )  {   # allow n outside CDS, 
					$cds_countT ++;					
					for my $k (0..length($seq)-1)  {
						   my $base = substr($seq, $k, 1);
						   $cdsT{$k}{$base} ++;
					  }
					print trueCDS ">$a[0]\_$a[1]\_$a[2]\_$a[3]\_$i\_CDS_$motif\_$max\_$frq\n$seq\n"; 

												
					if (exists $orf{"$a[0]\_$a[1]"} )   {
						my $codonseq = get_CDS_seq($a[2], $i, $a[3], \%genome);
						# get the sequence of codon before PAS
						
						my @pos = all_match_positions($codonseq,$orf{"$a[0]\_$a[1]"});
						next if scalar @pos != 1; # remove peaks that have multiple hits or no hit 
						my $cod_st = int($pos[0]/3) * 3;		# take the position of right codon, start within the selected range
						my $pas_DNA = substr($orf{"$a[0]\_$a[1]"}, $cod_st, ($aa_up + $aa_dn + 2)*3);

						#print "$motif\n$seq\n$codonseq\n$pas_DNA\n";
						my @pas_codon = CodonArray($pas_DNA, 3);
						print codonT "$motif\t$a[0]\_$a[1]\t$max\t$frq\t".join("\t", @pas_codon)."\n";  }
				 }
				else {	
					$intron_countT ++;
					for my $k (0..length($seq)-1)  {
						   my $base = substr($seq, $k, 1);
						   $cds_intronT{$k}{$base} ++;  }
					print trueIntron ">$a[0]\_$a[1]\_$a[2]\_$a[3]\_$i\_intron_$motif\_$max\_$frq\n$seq\n";  }   
			}  

			##### false signal count ######
			if (exists $false_pas{"$a[2]\t$i\t$a[3]"})   {

				my $motif = $false_pas{"$a[2]\t$i\t$a[3]"};
				my $seq = getSeq($a[2], $i, $a[3], \%genome);							
				# check whether the PAS region of peak is within CDS
				my $pas_st = $i;
				my $pas_end = $i + 6;
				if ($a[3] eq '-')  {
					$pas_st = $i-6;
					$pas_end = $i;}

				my $overlap = 0;
				for my $j ($pas_st .. $pas_end)  {
					$overlap ++ if !exists $cr{$j}; }
				
				if ( $overlap == 0 )  {   # allow n outside CDS, == 6 means all in coding region
					$cds_countF ++;					
					for my $k ( 0..length($seq)-1 )  {
						   my $base = substr($seq, $k, 1);
						   $cdsF{$k}{$base} ++;  }
					print falseCDS ">$a[0]\_$a[1]\_$a[2]\_$a[3]\_$i\_CDS_$motif\n$seq\n"; 

												
					if (exists $orf{"$a[0]\_$a[1]"} )   {
						my $codonseq = get_CDS_seq($a[2], $i, $a[3], \%genome);
						# get the sequence of codon before PAS
						my @pos = all_match_positions($seq,$orf{"$a[0]\_$a[1]"});
						next if scalar @pos != 1; # remove peaks that have multiple hits or no hit 
						my $cod_st = round($pos[0]/3, 0) * 3;		# take the position of right codon, start within the selected range
						my $pas_DNA = substr($orf{"$a[0]\_$a[1]"}, $cod_st, ($aa_up + $aa_dn + 2)*3);

						my @pas_codon = CodonArray($pas_DNA, 3);
						print codonF "$motif\t$a[0]\_$a[1]\t".join("\t", @pas_codon)."\n";  }

				 }
				else {	
					$intron_countF ++;
					for my $k (0..length($seq)-1)  {
						   my $base = substr($seq, $k, 1);
						   $cds_intronF{$k}{$base} ++;  }
					print falseIntron ">$a[0]\_$a[1]\_$a[2]\_$a[3]\_$i\_intron_$motif\n$seq\n";  }   				

			}   

	   	}
				   
	}

   ## end of search coding region 
   

   ## ========================================================================
   ## this part is for 3' UTR region, search 
   # define 3' UTR region
   my $st_utr; my $end_utr; 
	if ($a[3] eq '-')     {
		$st_utr = $a[4] - 100;
		$st_utr = 0 if $st_utr < 0;
		$end_utr = $a[6];}
	else   {
		$st_utr = $a[7];
		$end_utr = $a[5] + 100;}
	
	for my $i ($st_utr..$end_utr)   {
		##### false utr signal count ######
		if (exists $false_pas{"$a[2]\t$i\t$a[3]"})   {
			$utr_countF ++;
			my $motif = $false_pas{"$a[2]\t$i\t$a[3]"};

			my $seq = getSeq($a[2], $i, $a[3], \%genome);  
	
			for my $j (0..length($seq)-1)  {
				my $base = substr($seq, $j, 1);
				$utrF{$j}{$base} ++;  }
			print falseUTR ">$a[0]\_$a[1]\_$a[2]\_$a[3]\_$i\_UTR\_$motif\n$seq\n"; 
		 }


		#### true utr singal #####					     
		elsif (exists $true_pas{"$a[2]\t$i\t$a[3]"})   {
			$utr_countT ++;
			my($motif,$pos,$max,$frq) = split /\t/, $true_pas{"$a[2]\t$i\t$a[3]"};
			my $seq = getSeq($a[2], $i, $a[3], \%genome);  
			
			for my $j (0..length($seq)-1)  {
				my $base = substr($seq, $j, 1);
				$utrT{$j}{$base} ++;  }
			print trueUTR ">$a[0]\_$a[1]\_$a[2]\_$a[3]\_$i\_UTR\_$motif\t$max\t$frq\n$seq\n";  
		}
	}


	$item ++;
	if ($item == 500)  {
		$prog += $item;		
		$item = 0;
		print "$prog processed\n";
			}

	     } 
		
 
close REF;
close trueCDS;
close falseCDS;
close trueIntron;
close falseIntron;
close trueUTR;
close falseUTR;
close codonT;
close codonF;

#print Dumper(\%cdsT);
printGC(\%cdsT, $sam, "CDS_true");
printGC(\%cdsF, $sam, "CDS_false");
printGC(\%cds_intronT, $sam, "CDS_intron_true");
printGC(\%cds_intronF, $sam, "CDS_intron_false");
printGC(\%utrT, $sam, "UTR_true");
printGC(\%utrF, $sam, "UTR_false");


open (CT, ">$sam/classify/counts_TF_CDS_intron");
print CT "CDS_TRUE\tCDS_FALSE\tINTRON_TRUE\tINTRON_FALSE\tUTR_TRUE\tUTR_FALSE\n";
print CT "$cds_countT\t$cds_countF\t$intron_countT\t$intron_countF\t$utr_countT\t$utr_countF\n";

close CT;

## my %cdsT; my %cdsF; my %cds_intronT; my %cds_intronF;




##=============================================================================
## this part for subs
##=============================================================================



sub searchPeak     {
	my($chr, $start, $end, $ref) = @_;
	my %v; my @v; 
	for my $i ($start..$end)    {
		if (exists $$ref{$chr}{$i})   {
				my $val = $$ref{$chr}{$i};
				push @v, $val;
				$v{$val} = $i; }  }

	if (scalar @v > 0)  {
		my $max = max @v;
		my $sum = sum @v;
		my $pos = $v{$max}; 
		return ($pos,$max,$sum);
		     }
	else {
		return ($start, 0, 0); }
	}
				

				
sub getSeq     {
	my ($chr, $p, $strand, $ref) = @_;
	my $seq;
	if ($strand eq "-")   {
		$seq = substr($$ref{$chr}, $p-$dn, $dn+$length+6);
		$seq = reverse $seq;
		$seq =~ tr/ATCG/TAGC/; 
		}
	else   {
		$seq = substr($$ref{$chr}, $p-$length, $dn+$length+6);
		}
	return $seq;   }	
	

sub get_CDS_seq     {
	my ($chr, $p, $strand, $ref) = @_;
	my $seq;
	if ($strand eq "-")   {
		$seq = substr( $$ref{$chr}, $p-$aa_dn*3, ($aa_dn+$aa_up)*3+6 );
		$seq = reverse $seq;
		$seq =~ tr/ATCG/TAGC/; 
		}
	else   {
		$seq = substr($$ref{$chr}, $p-$aa_up*3, ($aa_dn+$aa_up)*3+6 );
		}
	return $seq;   }	
	



sub printGC   {

	my ($hash, $sample, $cat) = @_;
	my @nuc = ("A", "C", "G", "T");
	open (hand1, ">$sample/classify/gc_$cat.csv");
	print hand1 "pos\tA\tC\tG\tT\n";

	foreach my $i (sort keys %{$hash})  {
		my $p = $i-$length;
		print hand1 "$p";
		foreach my $j (@nuc) {
			if (exists $$hash{$i}{$j})  {
				print hand1 "\t$$hash{$i}{$j}";}
			else {
				print hand1 "\t0";}  }
		print hand1 "\n";  }
	
	close hand1;  }
	
sub CodonArray    {
	my ($seq, $len) = @_;
	my @array;
	for ( my $i = 0; $i <= length($seq)-$len; $i += 3 )  {
		my $cod = substr($seq, $i, $len);
		push @array, $cod; }
		return @array;    }
								



sub all_match_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /($regex)/g) {
        push @ret, pos($string)-length($1);
        #push @ret, [(pos($string)-length $1),pos($string)-1];  #this will report both start and end
    	}
 	return @ret; }



	
