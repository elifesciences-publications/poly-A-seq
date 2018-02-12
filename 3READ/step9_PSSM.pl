#!/usr/bin/perl -w
use List::Util qw(sum max);
#use List::MoreUtils qw(uniq);
use strict;
#use Getopt::Long;
use Data::Dumper;

my ($start, $end, $species);
$start = 0;
$end = 15; 
$species = "mouse_mm10";


mkdir "pssm";

# calculate the distribution of nucleotides in genome 
my $ge; my %genome; my %bkg_nuc;
open (hand2, "/home/yunkun/genomes/$species/Sequences/WholeGenomeFasta/genome.fa") or die $!;
while (<hand2>)   {
	$_=~ s/\s+$//;
	if (/^>/)       {
    		$ge=$_;
    		$ge=~ s/^>//;
    		next;}
	$_ =~ uc $_;
    	$genome{$ge} .= $_;
	$bkg_nuc{"A"} += $_ =~ tr/A//;
	$bkg_nuc{"C"} += $_ =~ tr/C//;
	$bkg_nuc{"T"} += $_ =~ tr/T//;
	$bkg_nuc{"G"} += $_ =~ tr/G//;
		}
close hand2;

my $total;
foreach (keys %bkg_nuc)  {
	$total += $bkg_nuc{$_}; }

foreach (keys %bkg_nuc)  {
	$bkg_nuc{$_} /= $total; }

print Dumper(\%bkg_nuc);

my @sample = ("SRR2225336", "SRR2225337", "SRR2225338", "SRR2225339",
           "SRR2225340", "SRR2225341", "SRR2225342", "SRR2225343", 
           "SRR2225344", "SRR2225345", "SRR2225346", "SRR2225347", 
           "SRR2225348", "SRR2225349", "SRR2225350", "SRR2225351");

foreach my $s (@sample)  {
	cal_score($s);   }
		
sub cal_score    {
	my ($sam) = @_;

	print "processing $sam\n";
	# the input sequence must be 101 nt long, with 50/51th nt are the cleavage		
	open (hand1, "$sam/polyA_UTR3.txt") or die $!;
	my @utr;
	while (<hand1>)  {
		chomp;
		next if /^>/;
		my $up = substr($_, (50+$start), $end - $start + 1 );
		push @utr, $up; }
	close hand1;	
	my %pssm = get_pwm(@utr);
	#print Dumper(\%pssm);

	# the sequences for Internal ORF and UTR
	
	#processing CDS files
	open (hand2, "$sam/polyA_CDS.txt") or die $!;
	open (out1, ">pssm/$sam\_CDS.txt");
	while (<hand2>)    {
		chomp;
		if (/^>/)  {
			$_ =~ s/>//;
			print out1 "$_\t";
			next }
		my $seq = substr($_, (50+$start), $end - $start + 1 );
		my $score = get_score($seq, \%pssm);
		print out1 "$score\n";  }
		
	close hand2; close out1;
	
	#processing UTR files 
	open (hand2, "$sam/polyA_UTR3.txt") or die $!;
	open (out1, ">pssm/$sam\_UTR3.txt");
	while (<hand2>)    {
		chomp;
		if (/^>/)  {
			$_ =~ s/>//;
			print out1 "$_\t";
			next }
		my $seq = substr($_, (50+$start), $end - $start + 1 );
		my $score = get_score($seq, \%pssm);
		print out1 "$score\n";  }
	close hand2; close out1;
		
		}

	



# ==========================================
# this part for supporting sub
# ==========================================


sub get_pwm {

    my @data = @_;
    my %pwm; my $line_number;
    foreach my $line (@data) {
	$line_number ++;
        for my $i (0 .. length($line) - 1)   {
		$pwm{substr($line, $i, 1)}{$i} ++;  
    }     }
	
	# get the weight of each entry
	foreach my $n (keys %pwm)   {
		foreach my $i ( keys %{$pwm{$n}} )  {
			if ( exists $bkg_nuc{$n} )  {
				$pwm{$n}{$i} = log2($pwm{$n}{$i}/($line_number*$bkg_nuc{$n}));   
			   }
			else {
				$pwm{$n}{$i} = "-Inf"; }
	}      }
 
    return %pwm;
}


sub get_score  {
	my ($seq,$ref) = @_;
	my $score = 0;
	for my $i (0..length($seq)-1)  {
		if (exists $$ref{substr($seq, $i, 1)}{$i} )  {
			$score += $$ref{substr($seq, $i, 1)}{$i};  }
		else {
			die {"no pssm entry was found!!"}   };  }
	return $score;
	}
	

sub log2 {
    my ($num) = @_;
    return log($num)/log(2);
    }
