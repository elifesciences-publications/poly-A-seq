#!/usr/bin/perl -w
#use List::Util qw(reduce max);
#use List::MoreUtils qw(uniq);
use strict;
#use Getopt::Long;


my %gene; my $ge;

open (hand1, "ref/CDS_DNA.fa") or die $!;

while (<hand1>)  {
	chomp;
	if (/^>/) {
		$ge = $_;
		$ge =~ s/>//;
		print "$ge\n";
		}
	else {
		$gene{$ge} .= $_;
    		}
}

close hand1; 

my @sample = ("SRR2225336", "SRR2225337", "SRR2225338", "SRR2225339",
           "SRR2225340", "SRR2225341", "SRR2225342", "SRR2225343", 
           "SRR2225344", "SRR2225345", "SRR2225346", "SRR2225347", 
           "SRR2225348", "SRR2225349", "SRR2225350", "SRR2225351");


foreach my $id (@sample)   {
	print "######### processing $id now ###########\n"; 
	count($id);  }


sub count   {

	my ($sam) = @_;


	my %can;
	open (hand1, "$sam/PAS_dicodons.txt") or die $!;

	while (<hand1>)     {
		next if /^Gene/;
		my @a = split /\t/;
		$can{$a[0]} += $a[2]; }
	close hand1;
		

	my %dicod; my %dicod_can; my %dicod_no;


	foreach my $g (keys %gene)    {
		for (my $i=0; $i<=length($gene{$g})-6; $i+=3)  {
			my $dicod = substr($gene{$g}, $i, 6);
			$dicod{$dicod} ++; 
			if (exists $can{$g})    {
				$dicod_can{$dicod} ++; }
			else {
				$dicod_no{$dicod} ++;  }
 }   }

	open (hand2, ">$sam/dicodon_genome_frequency.txt");
	print hand2 "dicodon\tgenome\tORF_cleavage\tnone\n";

	foreach my $di (keys %dicod)   {
		if (exists $dicod_can{$di} and exists $dicod_no{$di})  {
			print hand2 "$di\t$dicod{$di}\t$dicod_can{$di}\t$dicod_no{$di}\n";  }
		elsif (!exists $dicod_can{$di} and exists $dicod_no{$di} )  {
			print hand2 "$di\t$dicod{$di}\t0\t$dicod_no{$di}\n";  }
		elsif (exists $dicod_can{$di} and !exists $dicod_no{$di} )  {
			print hand2 "$di\t$dicod{$di}\t$dicod_can{$di}\t0\n";  }
		else {
			print hand2 "$di\t$dicod{$di}\t0\t0\n";  }    }
	close hand2;    }




