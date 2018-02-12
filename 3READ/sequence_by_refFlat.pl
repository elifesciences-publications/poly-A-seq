#！usr/bin/perl -w
use List::Compare;
use strict;
#use Bio::Util::codonUsage qw(translate);

=head1 the purpose
This script is to extract the sequence from refFlat file, which can be created with gtf or gff with UCSC tools gtfToGenePred
It is still unknown why some sequence give wrong coding region, but these sequences, or id will be stored in a translate_warning files, with original DNA sequence and translation. 
As a refrence, this warning file will also be taken to excluding them in the further analyses of PRF, which in this PolyA seq case, are CDS signal.
For codon frequency, these warning data are not included as well.  

=cut

my $species = "NC10";

mkdir "ref";
my $ge; my %genome;
open (hand2, "/home/yunkun/genomes/$species/Sequences/WholeGenomeFasta/genome.fa") or die $!;
while (<hand2>)   {
	$_=~ s/\s+$//;
    if (/^>/)       {
    $ge=$_;
    $ge=~ s/^>//;
    next;}
    $genome{$ge} .= $_; }

close hand2;
	
# use refflat reference file the annotate the peaks, either at 3' end or internal CDS	
my %cds; my %utr; my %cod;
open (REF, "/home/yunkun/genomes/$species/Genes/refFlat.txt") or die $!;
open (Seq, ">ref/CDS_DNA.fa");
open (P, ">ref/CDS_pep.fa");
open (W, ">ref/translate_warning.txt");
open (Size, ">ref/CDS_intron_size.txt");
print Size "gene_id\tlocus\ttranscript\tCDS\tCDS_intron\n";
open (cds_exon, ">ref/CDS_exons.txt");

while (<REF>)   {
	
   chomp;
   my @a = split /\t/;
   my @ins = split /,/, $a[9];
   my @ine = split /,/, $a[10];
   my $cds_count = 0;
   my $utr_count = 0;
   print "processing $a[0]\n";
   # define CDS region, first get all exon
   my @cds;
    if ($a[7] - $a[6] > 0)   {
	   
		for my $i ( 0..$a[8]-1 )  {
			
			for my $j ($ins[$i]..$ine[$i]-1)  {
				
				push @cds, $j; }  }
		my $full = @cds;
		#print "all=".join(",",@cds)."\n";
		my @utrL; my @utrR;
		@utrL = $a[4]..$a[6]-1 if $a[4] - $a[6] != 0;
		@utrR = $a[7]..$a[5] if $a[5] - $a[7] != 0;
		my $lcl = List::Compare->new(\@cds, \@utrL);
		@cds = $lcl->get_unique;
		# needs to sort by value, not sort by 
		@cds = sort {$a <=> $b} @cds;
		#print join(",",@cds)."\n";
		my $lcr = List::Compare->new(\@cds, \@utrR);
		@cds = $lcr->get_unique;
		@cds = sort {$a <=> $b} @cds;
		#print join(",",@cds)."\n";
		my $CDS = @cds;
		my $cds_intron = $a[7]-$a[6]-$CDS;
		print Size "$a[0]\t$a[1]\t$full\t$CDS\t$cds_intron\n";
			
	my $seq; 
	my %cds_exon;
	my $position = $cds[0];
	my $ex_number = 1;

	foreach my $i (@cds)  {
		if ($i <= $position + 1)  {
			#print "$i,$position\t";
			$cds_exon{$ex_number} .= substr($genome{$a[2]}, $i, 1) ;
			$position ++;  }
		else {
			$ex_number ++;
			#print "\n\n";
			$cds_exon{$ex_number} .= substr($genome{$a[2]}, $i, 1) ;
			$position = $i;  }
		
		$seq .= substr($genome{$a[2]}, $i, 1) ;
		$seq =~ tr/atcg/ATCG/;  #some sequence are in lower case
		}
		
	if ($a[3] eq '-')  {
		foreach my $id (keys %cds_exon)  {
			$cds_exon{$id} = reverse $cds_exon{$id};
			$cds_exon{$id} =~ tr/ATCG/TAGC/;  }

		$seq = reverse $seq;
		$seq =~ tr/ATCG/TAGC/;  }

	#report the cds exons
	
	if ( $a[3] eq "-")  {
		my $exon_order = 1;
		foreach my $id (sort {$b <=> $a} keys %cds_exon)  {
			print cds_exon ">$a[0]\_$a[1]\_exon_$exon_order\n$cds_exon{$id}\n";
			$exon_order ++; } 
	}
	else {
		foreach my $id (sort {$a <=> $b} keys %cds_exon)  {
			print cds_exon ">$a[0]\_$a[1]\_exon_$id\n$cds_exon{$id}\n"; }
	}


	my $protein; 
	for ( my $i = 0; $i < length($seq)-1; $i += 3 )    {
		my $aa = translate(substr($seq, $i, 3));
		$protein .= $aa; }

	if ($protein !~ /X|^.+\*.+$/ and $protein =~ /^M.+\*$/)  {
	
		print Seq ">$a[0]\_$a[1]\n$seq\n";
		print P ">$a[0]\_$a[1]\n$protein\n";
		# count codons in proteins that are 
		# 1, started with ATG and stop with stop codons
		# 2, no internal * and no irrugular X, here X are not right codon, either cotain N or not 3 nt
		for ( my $i = 0; $i < length($seq)-1; $i += 3 )    {
			my $cod = substr($seq, $i, 3);
			$cod{$cod} ++; }      }
	else {	
		print W ">$a[0]\_$a[1]\n$protein\n$seq\n";}
	}  }
		
close REF; close Seq; close cds_exon;	

open (hand4, ">ref/codon_frequency.txt");

foreach my $c (keys %cod)   {
	my $aa = translate($c);
	print hand4 "$c\t$aa\t$cod{$c}\n";}

close hand4;














sub translate {
 my ($cod) = shift;
 
 my (%codon2aa) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );
 
 if (exists $codon2aa{$cod}) {
  	return $codon2aa{$cod};  } 
 else { 
	print "warning, $cod is not a codon\n";
	return "X";}
              }




