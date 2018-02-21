#！usr/bin/perl -w
use List::MoreUtils qw(uniq);
use strict;

# Author: Yunkun Dang
# Updated on Feb 18, 2018


my $ch; my %genome;

### put the path and name of your genome fasta data. must including all chromosomes in one file. 
open (hand1, "~/genomes/mouse_mm10/Sequence/WholeGenomeFasta/genome.fa") or die $!;


while (<hand1>)   {
	$_=~ s/\s+$//;
    if (/^>/)       {
    $ch=$_;
    $ch=~ s/^>//;
    next;}
    $genome{$ch} .= $_; }

close hand1;

### put your sample name here, quote every name and seperate with comma. sample name are the same as used in step1
my @sample = ("SRR2225336", "SRR2225337", "SRR2225338", "SRR2225339", "SRR2225340", "SRR2225341", "SRR2225342", "SRR2225343", "SRR2225344", "SRR2225345", "SRR2225346", "SRR2225347", "SRR2225348", "SRR2225349", "SRR2225350", "SRR2225351");

foreach my $id (@sample)   {
	print "processing $id now\n\n";
	filter($id);  }


sub filter       {
	my($sam) = @_;
	open (hand2, "$sam/thout/$sam\_uniq.bed") or die $!;
	open (posout, ">$sam/thout/$sam\_filtered_pos.bed");
	open (out, ">$sam/thout/$sam\_filtered.bed");
		
	while (<hand2>)     {
		chomp;
		my @a = split /\t/;
		my @b = split /_A=/, $a[3];
		# retrive the sequence from genome for region after the 3' end
		my $seq;
		if ($a[5] eq '-')   {
			$seq = substr($genome{$a[0]}, $a[1]-$b[1], $b[1]);
			$seq = reverse $seq;
			$seq =~ tr/ATGC/TACG/;  }
		else {
			$seq = substr($genome{$a[0]}, $a[2], $b[1]); }
		
		if ($seq =~ /[TCG]{2}/)    {
			print out "$_\n";
			#print "$seq\n";
			if ($a[5] eq '-')  {
				my $end = $a[1]+1;
				print posout "$a[0]\t$a[1]\t$end\t$a[3]\t$a[4]\t$a[5]\n";}
			else   {
				my $start = $a[2] - 1;
				print posout "$a[0]\t$start\t$a[2]\t$a[3]\t$a[4]\t$a[5]\n";}
			}	
	
		
		}
		
	close hand2;  }
			
