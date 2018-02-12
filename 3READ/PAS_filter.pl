#！usr/bin/perl -w
use strict;

=head1 The Methodology to anthenticate the real PAS reads
Author: Yunkun 
date: May 5, 2017

Only apply for this specific case. i.e. primed with XXXT12VN primer (2P seq) and sequence 
the reads from 5' end (universal primer)
The way to do this is 
1, pile up the 3' end (putative cutting site) with bedtools
2, look at each position with 2 filter criteria: 
	1, support by < 3 read 
	2, have >= 6 consecutive A or 7 A within -10 to +10 (a 20nt window around cut site)
	
=cut

#my $species = "mouse_mm10";


use Getopt::Long;


my ($sam, $strand, $species);

$sam = "";
$strand = "plus";
$species = "NC10";

GetOptions 
(
  's=s'=>\$sam,		# name of sample
  'strand=s'=>\$strand, 	# strand
  'species=s'=>\$species	# species
  );

print "============processing $sam =============\n";

my $ch; my %genome;
open (hand1, "/home/yunkun/genomes/$species/Sequences/WholeGenomeFasta/genome.fa") or die $!;
while (<hand1>)   {
	$_=~ s/\s+$//;
    if (/^>/)       {
    $ch=$_;
    $ch=~ s/^>//;
    next;}
    $genome{$ch} .= $_; }

close hand1;

open (hand2, "$sam/thout/pos_$strand\_p.bedgraph") or die $!;
open (out, ">$sam/thout/pos_$strand.bedgraph");

my %pos;  #record the position that are filtered, i.e. the false signal
my $count = 0; # count number of false reads	
my $total = 0; # count total number 	
while (<hand2>)     {
	chomp;
	my @a = split /\t/;
	#if ($a[3] >= 5)   {
		#print out "$_\n";
		#next } 
	$total += $a[3];	
	#get sequence around the putative cleavage
	
	my $seq;	
	# deal with two criteria
	
	if ($strand eq "neg")   {
		$seq = substr($genome{$a[0]}, $a[1]-20, 20);
		$seq =~ tr/atcg/ATCG/;
		$seq = reverse $seq;
		$seq =~ tr/ATCG/TAGC/}
	else {
		$seq = substr($genome{$a[0]}, $a[1], 20);}
		$seq =~ tr/atcg/ATCG/; 

	#count the number of A
	#my $d = () = ($n =~ /A/g);
	
	my $on = 0; # the switch

	Search: for my $i (0..length($seq)-12)  {
		my $fragment = substr($seq, $i, 12);
		my $d = ($fragment =~ tr/A//); # this is faster
		if ($seq =~ m/AAAAAA|A{3,}.A{2,}.A{2,}/g or $d >= 8) {
			$pos{"$a[0]\_$a[1]"} = ""; 
			$count += $a[3];
			$on = 1;
			last Search; }   }

	print  out "$_\n" if $on == 0; 
	
	  }
	
close hand2; close out;	

# filter the bed file
open (hand3, "$sam/thout/$sam\_filtered_pos.bed") or die $!;
open (out1, ">$sam/thout/$strand\_filtered_PAS.bed");

my $n = 0; 
my $std = "+"; 
$std = "-" if $strand ne "plus";

while (<hand3>)    {
	chomp;
	my @a = split /\t/;
	
	next if $a[5] ne $std;

	if (exists $pos{"$a[0]\_$a[1]"})  {
		$n ++; 
		next }

	print out1 "$_\n";  }


close hand3; close out1;
print "filtered $std strand reads (total = $total) from $sam = $count; from bed file is $n\n";














		
