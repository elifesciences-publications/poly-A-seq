#!/usr/bin/perl -w
#use List::Util qw(reduce max);
#use List::MoreUtils qw(uniq);
use strict;
use Getopt::Long;
#use List::Compare;
#use Bio::Util::codonUsage qw(translate);

=head1
First created: August 10, 2017
Author = Yunkun Dang
purpose: 
looking for peaks (i.e. cleavage site) that are locationed inside the coding regions (ORF). Find out the upstream region (default is 30 to 10 upstream of cleavage site), report the codons inside in these regions. Note that the selected upstream region must be inside the coding region. Meanwhile, also randomly pick a site within selected genes (i.e. the same genes that have the internal cleavage sites) and report the sequences of codons. 


=cut




my ($sam, $genome, $refflat, $upstream, $downstream);

$sam = "";
$genome = "";
$refflat = "";
$upstream = 30;
$downstream = 10;


GetOptions 
(
  's=s'=>\$sam,		# name of sample
  'g=s'=>\$genome,	# path and file name of genome
  'ref=s'=>\$refflat,		# refflat file, can be generated from gtf or gff 
  'start=s'=>\$upstream,	# the upstream limit for search the region of PAS
  'end=s'=>\$downstream, 	# downstream limit for PAS
  );

print "============processing $sam =============\n";

if ($sam eq "" or $genome eq "" or $refflat eq "") {
	die "Missing one or more input file name(s)!
	use of script:
	perl findPAS.pl -s <sample_name> -g <path to genome fasta file ) -ref <refflat file> -start <num> -end <num>\n"} # check parameters

# loading reference of cleavage site with density information
my %plus; my %minus;
open (hand1, "$sam/thout/pos_plus.bedgraph") or die "$sam pos_plus.bedgraph not in thout/";
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

# refernce genome sequences and ORF sequence, exlude the genes that are labeled as ORF but have internal stop codons
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

close hand2;  close hand3;

my %cai;
open (hand3, "ref/CBI_ref.txt") or die "no reference CAI/CBI file";
while (<hand3>)   {
	$_=~ s/\s+$//;
	next if /^AA/;
	my ($aa,$codon,$cai,$rand) = split /\t/;
    $cai{$codon} = $cai; }
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
	
	## IMPORTANT: avoid using hash in hash unless it is very necessary, because it would cause circular reference and comsume huge amount of RAM.	
	$peak{"$chr\t$p\t$strand"} = "$max\t$frq";
	}
	
# use refflat reference file the annotate the peaks, either at 3' end or internal CDS	
my %cds; my %ran_cds; my %diCod; my %ran_diCod;  

open (REF, $refflat) or die $!;
open (codon, ">$sam/PAS_codons.txt");
open (codFA, ">$sam/PAS_codons.fa");
open (dicodon, ">$sam/PAS_dicodons.txt");
open (ran_cod, ">$sam/PAS_random_codons.txt");
open (ran_dicod, ">$sam/PAS_random_dicodons.txt");
open (CAI, ">$sam/PAS_CAI.txt");
open (ran_CAI, ">$sam/PAS_random_CAI.txt");
open (AA, ">$sam/PAS_amino_acid.fa");

my $item = 0;
my $prog = 0;

while (<REF>)   {
   chomp;
   my @a = split /\t/;
   my @ins = split /,/, $a[9];
   my @ine = split /,/, $a[10];


   # define coding region, first get all exon
   my %cr;
    if ($a[7] - $a[6] > 0 and !exists $exc{"$a[0]\_$a[1]"} )   {
	   
		for my $i ( 0..$a[8]-1 )  {
			for my $j ($ins[$i]..$ine[$i]-1)  {
				$cr{$j} = ""; }  }
		
		# remove utrL from hash cr
		if ($a[4] - $a[6] != 0)   {
			for my $i ($a[4]..$a[6]-1)  {
				delete $cr{$i} if exists $cr{$i}; } }
		
		# remove utrR from hash cr		
		if ($a[5] - $a[7] != 0)   {
			for my $i ($a[4]..$a[6]-1)  {
				delete $cr{$i} if exists $cr{$i}; } }
		

		# searching the coding region, including intron, and found whether the 
		# PAS region is within the coding region (i.e, the -30 to -10 is completely
		# inside coding region.
		foreach my $i ($a[6]..$a[7])  {
			if (exists $peak{"$a[2]\t$i\t$a[3]"})   {
				my($max,$frq) = split /\t/, $peak{"$a[2]\t$i\t$a[3]"};

				# remove peaks with low reads
				#next if $frq < 5;				
				# check whether the PAS region of peak is within CDS
				my $pas_st = $i-$upstream;
				my $pas_end =$i-$downstream;
				if ($a[3] eq '-')  {
					$pas_st = $i+$downstream;
					$pas_end = $i+$upstream;}
				
				# search if desired PAS region is inside
				my $overlap = 0;
				for my $j ($pas_st .. $pas_end)  {
					$overlap ++ if !exists $cr{$j}; }
				
				if ( $overlap == 0 )  {   # allow n outside CDS
					
						my $seq = getSeq($a[2], $pas_st, $pas_end, $a[3], \%genome);
												
						if (exists $orf{"$a[0]\_$a[1]"} )   {

							# get the random codon sequence (in-frame) with desired length, take 10 times each 
							my $ran_pos; my $ran_seq;
							for (0..9)   {
								$ran_pos = abs int(rand(length($orf{"$a[0]\_$a[1]"})/3)*3 - int(($upstream-$downstream)/3+1)*3);
								$ran_seq .= substr($orf{"$a[0]\_$a[1]"}, $ran_pos, int(($upstream-$downstream)/3)*3);  }
							
							my @ran_cai = CodonArrayCAI($ran_seq, \%cai); 
							my @ran_codon = CodonArray($ran_seq, 3);
							my @ran_diCod = CodonArray($ran_seq, 6);
							
							foreach my $id (@ran_diCod)  {
								$ran_diCod{$id} ++; }
								
							foreach my $id (@ran_codon)  {
								$ran_cds{$id} ++; }
							
							my $dis = ($upstream-$downstream)/3;
							for my $n (0..9)    {
								my @ran_codon_sub = splice(@ran_codon, 0, $dis);
								my @ran_cai_sub = splice(@ran_cai, 0, $dis);
								my @ran_diCod = splice(@ran_diCod, 0, $dis-1);	
								print ran_cod "$a[0]\_$a[1]\t".join("\t", @ran_codon_sub)."\n";
								print ran_CAI "$a[0]\_$a[1]\t".join("\t", @ran_cai_sub)."\n"; 
								print ran_dicod "$a[0]\_$a[1]\t".join("\t", @ran_diCod)."\n"; }
							
							# get the sequence of codon before PAS
							my @pos=all_match_positions($seq,$orf{"$a[0]\_$a[1]"});
							next if scalar @pos > 1; # remove peaks that have multiple hits   
							my $cod_st = int($pos[0]/3 + 1) * 3;		# take the position of right codon, start within the selected range
							my $pas_DNA = substr($orf{"$a[0]\_$a[1]"}, $cod_st, int(($upstream-$downstream)/3)*3);

							print codFA ">$a[0]\_$a[1]\_$i\n$pas_DNA\n";
							
							my @pas_cai = CodonArrayCAI($pas_DNA, \%cai); 
							my @pas_codon = CodonArray($pas_DNA, 3);
							my @pas_diCod = CodonArray($pas_DNA, 6);

											
							foreach my $id (@pas_diCod)  {
								$diCod{$id} ++; }
								
							my $pas_pep;	
							foreach my $id (@pas_codon)  {
								$cds{$id} ++;
								$pas_pep .= translate($id); }
							
				
							print codon "$a[0]\_$a[1]\t$max\t$frq\t".join("\t", @pas_codon)."\n";
							print CAI "$a[0]\_$a[1]\t$max\t$frq\t".join("\t", @pas_cai)."\n";
							print AA ">$a[0]\_$a[1]\_$i\n$pas_pep\n";
							print dicodon "$a[0]\_$a[1]\t$max\t$frq\t".join("\t", @pas_diCod)."\n";
						}

					}   }
		
		}      }
				   
	
   ##--------------------------------------------------------------------------------------
   # show progress
	$item ++;
	if ($item == 500)  {
		$prog += $item;		
		$item = 0;
		print "$prog processed\n";
			}

	     } 


close REF; close codon; close dicodon; close ran_cod; close ran_dicod;
 close CAI; close AA; close ran_CAI; 

##==============================================================
## this part for subs

sub CodonArray    {
	my ($seq, $len) = @_;
	my @array;
	for ( my $i = 0; $i <= length($seq)-$len; $i += 3 )  {
		my $cod = substr($seq, $i, $len);
		push @array, $cod; }
		return @array;    }
								
								
sub CodonArrayCAI    {
	my ($seq, $ref) = @_;
	my @array;
	for ( my $i = 0; $i < length($seq); $i += 3 )  {
		my $cod = substr($seq, $i, 3);
		push @array, $$ref{$cod}; }
		return @array;    }


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
	my ($chr, $st, $e, $strand, $ref) = @_;
	my $seq = substr($$ref{$chr}, $st, $e-$st);
	if ($strand eq "-")   {
		$seq = reverse $seq;
		$seq =~ tr/ATCG/TAGC/; }
	return $seq;   }	
	


sub match_positions {
    my ($regex, $string) = @_;
    return if not $string =~ /($regex)/;
    return (pos($string), pos($string) + length $1);
}


sub all_match_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /($regex)/g) {
        push @ret, pos($string)-length($1);
        #push @ret, [(pos($string)-length $1),pos($string)-1];  #this will report both start and end
    }
    return @ret
}
