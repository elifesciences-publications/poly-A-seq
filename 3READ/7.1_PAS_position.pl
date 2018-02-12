#!/usr/bin/perl -w
use strict;
use Getopt::Long;

=head1
this script is meant search the position of PAS signal based on the criteria
The pattern used for the search will be from a file generated in step5

=cut

my ($sam, $genome, $reference, $num);

$sam = "";
$genome = "";
$reference = "";
$num = 5; # default number to select PAS singal sequences for analyses 



GetOptions 
(
  's=s'=>\$sam,		# name of sample
  'g=s'=>\$genome,	# path and file name of genome
  'ref=s'=>\$reference,	# reference file to store the found PAS by step6
  'tr=s'=>\$num		# number of PAS motifs for search position
  );

if ($sam eq "" or $genome eq "" or $reference eq "") {
	die "Missing one or more input file name(s)!
	use of script:
	perl findPAS.pl -s <sample_name> -g <path to genome fasta file ) -ref <PAS reference file> -tr <number of PAS sequence>\n"} # check parameters


my $ge; my %genome;
open (hand2, $genome) or die $!;
while (<hand2>)   {
	$_=~ s/\s+$//;
	if (/^>/)       {
    		$ge=$_;
    		$ge=~ s/^>//;
    		next;}
	$_ =~ uc $_;
    	$genome{$ge} .= $_; }

my $pattern; 

# motif files from analyses
open (hand1, "$reference") or die $!;

# output file, to a bed format
open (hand3, ">$sam/classify/polyA_prediction_plus.bed");
open (hand4, ">$sam/classify/polyA_prediction_minus.bed");

while (<hand1>)  {
	next if /motif/;
	last if $num == 0;
	my @a = split /\t/;
	$pattern .= "$a[0]\|";
	# we assume that it would have a T rich region before
	$num --; }

$pattern =~ s/\|$//;
$pattern =~ s/U/T/g;
print $pattern."\n";
# define the pattern


foreach my $i (sort keys %genome)  {
	my @pos_plus = all_match_positions($pattern, $genome{$i});
	my $genome_bottom = reverse($genome{$i});
	$genome_bottom =~ tr/ATCG/TAGC/;
	my @pos_minus = all_match_positions($pattern, $genome_bottom);


	my $n = 0; my $e = 0;	
	foreach my $p (@pos_plus) {
		print hand3 "$i\t$$p[0]\t$$p[1]\t$$p[2]\_+\t.\t+\n";
			}
	

	foreach my $p (@pos_minus) {
		my $end = length($genome{$i}) - $$p[0];
		my $begin = $end - length($$p[2]);
		print hand4 "$i\t$begin\t$end\t$$p[2]\_-\t.\t-\n";
			}

		     }



close hand2; close hand3; close hand4;

## ================================================================================================
## This section for sub ##
## ================================================================================================
sub match_positions {
    my ($regex, $string) = @_;
    	if ($string =~ /($regex)?/)   {
    		return (pos($string) - length $1, pos($string)); }
	else {
		return}
		
}

sub all_match_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /($regex)/g) {
        push @ret, [(pos($string)-length $1), pos($string) ,$1];
    }
    return @ret
}




