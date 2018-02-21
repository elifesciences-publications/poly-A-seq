#!usr/bin/perl -w
use strict;

=head1 The Methodology to anthenticate the real PAS reads
Author: Yunkun 
date: Feb 15, 2018

Based on the previous analyses, we know each reads has 1-n A tailing the read,
so after the mapping we are looking back the A tail and check whether these tails 
are really untemplate A or just simply came from sequnces contain A at the 3' end

=cut

#my $species = "NC10";

my $ch; my %genome;

# reference genomic sequence data here.
open (hand1, "/home/yunkun/Deepseq/Indexes/NC10.fa") or die $!;


while (<hand1>)   {
	$_=~ s/\s+$//;
    if (/^>/)       {
    $ch=$_;
    $ch=~ s/^>//;
    next;}
    $genome{$ch} .= $_; }

close hand1;

# sample name here. should be the same as used in step1
my @sample = ("CHX_1", "CHX_2", "Nuc_1", "Nuc_2", "Total_1", "Total_2");

foreach my $id (@sample)   {
	print "######### processing $id now ###########\n"; 
	filter($id);  }

sub filter       {
	my($sam) = @_;
	open (hand2, "$sam/thout/$sam\_uniq.bed") or die $!;
	open (posout, ">$sam/thout/$sam\_filtered_pos.bed");
	open (out, ">$sam/thout/$sam\_filtered.bed");
		
	while (<hand2>)     {
		chomp;
		my @a = split /\t/;
		# deal with two conditions, this one is for 3P seq
		if ($a[3] =~ /_A=/)   {
			my @b = split /_A=/, $a[3];
			# retrive the sequence from genome for region after the 3' end
			my $seq;
			if ($a[5] eq '-')   {
				$seq = substr($genome{$a[0]}, $a[1]-$b[1], $b[1]);
				$seq =~ tr/atcg/ATCG/;
				$seq = reverse $seq;
				$seq =~ tr/ATGC/TACG/;  }
			else {
				$seq = substr($genome{$a[0]}, $a[2], $b[1]);
				$seq =~ tr/atcg/ATCG/; }
			# check untemplate A, at least one nucleotide other than A
			if ($seq =~ /[TCG]{1,}/)    {
				print out "$_\n";
				my $anno = $a[3]."_$seq";
				#print "$seq\n";
				if ($a[5] eq '-')  {
					my $end = $a[1]+1;
					print posout "$a[0]\t$a[1]\t$end\t$anno\t$a[4]\t$a[5]\n";}
				else   {
					my $start = $a[2] - 1;
					print posout "$a[0]\t$start\t$a[2]\t$anno\t$a[4]\t$a[5]\n";}
			}       }
		# THIS ONE IS FOR 2P-SEQ, no untemplate A at 3'end
		else {
			if ($a[5] eq '-')  {
				my $end = $a[1]+1;
				print out "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t-\n";
				print posout "$a[0]\t$a[1]\t$end\t$a[3]\t$a[4]\t$a[5]\n";}
			else   {
				my $start = $a[2] - 1;
				print out "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t+\n";
				print posout "$a[0]\t$start\t$a[2]\t$a[3]\t$a[4]\t$a[5]\n";}
			}       
	
	
		
		}
		
	close hand2;  }





			
