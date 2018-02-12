#！usr/bin/perl -w

use List::Util qw(shuffle max reduce);
use strict;
use List::MoreUtils qw(uniq);
#use Statistics::R;
use Math::CDF qw(pbinom);
#use Number::Format qw(round);
use Getopt::Long;
#use Cwd 'abs_path';
use Data::Dumper;

my ($start,$end,$file_in,$file_out,$Pval,$R);

$start = -30; # default values
$end = -10;
$file_in="";
$Pval = 0.05;
$R = 0.01;


GetOptions 
(
  'start=s'=>\$start,
  'end=s'=>\$end,
  'in=s'=>\$file_in,
  'P=s'=>\$Pval,    			# threshold Pval of binomial distribution
  'ratio=s'=>\$R				# ratio of hits to total sequences
);

if ($file_in eq "") {
	die "Missing input file name(s)!
	use of script:
	perl findPAS.pl -start -30 -end -10 -in <file_out> -P 0.05 -ratio 0.01\n"} # check parameters

#$file_in =~ s/^.+\///g;
my $file_out = $file_in;
$file_out =~ s/\.[a-z]+//;

my %can;	# candicate sequences
my $can;
my %clv;	# cleavage preference

open (hand1, $file_in) or die "no input file found";


while (<hand1>)     {
	chomp;
	if (/^>/)   {
		$can = $_;
		$can =~ s/^>//;
		next }
	my $seq = substr($_, int(length($_)/2) + $start, abs($start - $end));
	my $c = substr($_, int(length($_)/2), 2);
	$clv{$c}++;
	$can{$can} = $seq;
}	

close hand1;

print Dumper(\%clv);

# get random sequence from the dinucleotide array. and substring 20nuc randomly from the sequence
my $ratio = 1.0;

open (hand2, ">$file_out\_$start-$end\_sig_motifs.txt");
print hand2 "motif\tfrqency\tfrequency_peak\ttotal\trand_frq\trand_frq_peak\tPval_frq\tPval_frq_peak\n";
 
Search:	while ($ratio > 0.01)   {
			
			my %hex =();		#number of hexamer 
			my %hexNum;			# number of sequence contain that pattern (>=1)
			my @ranDi = ();		#random dimer array
			my $k = scalar keys %can; #number used in binomial ditribution as number of test
			# get the hexamer and dimer
			foreach my $id (keys %can)    {
				
				# count hexamer
				for my $i (0..length($can{$id})-1-6)   {
					my $m = substr($can{$id}, $i, 6);   #print "$m\n";
					$hex{$m} ++; 
					push @{$hexNum{$m}}, $id; }  			#record id, which reflect the number of sequence contain the hexamer
					
				# count dimer
				for my $i (0..length($can{$id})-1-2)   {
					my $m = substr($can{$id}, $i, 2);   #print "$m\n";
					push @ranDi, $m;  }  
				}	
			
			### set Pval < 0.05 for the highest frequency of hexamer, by comparing the random sequences that containing the hexamer
			my $p1 = 1;
			
	Pval:	while ($p1 > $Pval)    {
				# get the key of element of %hex with highest value     
				my $hex_highest = List::Util::reduce { $hex{$b} > $hex{$a} ? $b : $a } keys %hex;
				my @n_seq = uniq(@{$hexNum{$hex_highest}}); 	# the number of sequences that contain >= 1 times of hexamer
				my $n_seq = scalar @n_seq;
								
				# terminate Search loop upon < 1% sequence contain that hexamer
				$ratio = $n_seq/$k;
				last Search if $ratio < 0.01;	
				
				my ($ranCount, $ranFrq) = ranFrq($hex_highest, \@ranDi);

				my $fold = 0;				
				if ( $ranFrq != 0 )  {
					$fold = $n_seq/$ranFrq;  }
 
				my $p2 = 1 - pbinom($hex{$hex_highest}-1, $k, $ranCount/1000);
				my $p1 = 1 - pbinom($n_seq-1, $k, $ranFrq/1000);
				#print hand2 "$hex_highest\t$hex{$hex_highest}\t$n_seq\t$k\t$ranCount\t$ranFrq\t$p1\t$p2\n"; 

				if ($p1 < 0.05 && $fold >= 2 )    {
					my $Rmotif = $hex_highest;
					$Rmotif =~ tr/T/U/;
					print hand2 "$Rmotif\t$hex{$hex_highest}\t$n_seq\t$k\t$ranCount\t$ranFrq\t$p1\t$p2\n"; 
					# remove the keys from %can that contain $hex_highest
					foreach (keys %can)    {
						delete $can{$_} if $can{$_} =~ m/$hex_highest/;   }
						last Pval;  }
				else {
					delete $hex{$hex_highest};   }
				}
			}

	



close hand2; 
	    



 


##=========================================================================
## This part is for subs
##=========================================================================

# function 1, count the number of pattern in 1000 random sequences that made from  
sub ranFrq      {
	my($pattern, $diSeq) = @_;

	# get random sequence from the dinucleotide array. by randomly pick elements, put in an array called @sanSeq
	my @ranSeq; 	#random sequences
	
	for (0..999)   {
		my $ranSeq;
		for (0..9)   {
			my $ran = abs(int(rand(scalar @$diSeq))-1);
			$ranSeq .= @$diSeq[$ran];  }
		push @ranSeq, $ranSeq; }
	
	# count the number of matched hexamer in @ranSeq, freqency means number of sequence that at least have one pattern
	my $ranCount;	my $frq;
	foreach my $seq (@ranSeq)   {	
		my @matches = ($seq =~ /$pattern/g);
		my $d = @matches;
		$frq ++ if $d > 0;			
		$ranCount += $d;   }
	
	return ($ranCount, $frq);  }
		

__END__

This part could be used later to connect R, but it runs quite slow.

			# Pass and retrieve data (scalars or arrays) with Statistics::R
				# use the R function binomial distribution
				my $R = Statistics::R->new();
				$R->set('n', $n_seq);
				$R->set('k', $k);
				$R->set('p', $ranFrq/1000);
				$R->run(q`y <- 1-pbinom(n-1,k,p)`);
				my $p = $R->get('y');
				print "pVal = $p\n";
				$R->stop();



