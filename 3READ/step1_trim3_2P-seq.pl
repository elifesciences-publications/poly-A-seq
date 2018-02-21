use strict;
use warnings;

# Feb 5, 2018
# Yunkun Dang


# set minimum number of A
my $min=10;
# set minimum length of reads after trimming poly(A) stretch
my $size = 20;

=head1
IMPORTANT: 
These sample were done with 2P-seq. However, due to the mistake of using forward primer as sequencing primer, the result shows the sequence as XXXXAAAAAAAA...., so the method for trimming will be very different with other 2P seq protocol for analyses. 
Put your fastq files into a folder named "fq" 
=cut


# sample name here. THe file name will be, for example, Nuc_1.fastq)
my @sample = ("Nuc_1", "Nuc_2", "Total_1", "Total_2");


foreach my $id (@sample)   {
	print "######### processing $id now ###########\n"; 
	trim3($id);  }



####trim 3'
sub trim3                                                      {
my ($sam)=@_; my $count=0; my $id=10000000; my $idfq; my $on=0; my $si=0; my $seq;
my $cca=0;my %cca=(); my $t = 0;

mkdir $sam;
mkdir 'fq';
open (hand1,"fq/$sam.fastq") or die $!;
open (hand2,">$sam/$sam.txt");
open (hand3,">fq/$sam\_trimA.fq");

while (<hand1>)                 {
	$_ =~ s/\s+$//;
	$count++;
	my $mod=$count % 4;
	my $len = length($_);

	if ($mod == 1)       {
		$idfq=$_; 
		$idfq =~ s/ /_/g;
		next  }

	if ($mod == 2)             {
		$on=0; $si=0;
	if (/^(.+?)A{$min,}/g)   {
		$seq=$1; #print "$seq\n";
		 }
	else {next}
	

	if (length($seq) >= $size)   {
		$on=1; $id++;
		$si=length($seq);
		print hand2 $seq,"\n";
		print hand3 "$idfq\n";
		print hand3 "$seq\n";
		print hand3 "+$id\n";}     }
	elsif ($mod ==0 && $on==1) {
		my $qua=substr($_,0,$si);
		print hand3 $qua,"\n";     }   }
close hand1; close hand2; close hand3; 


#### uni
open (hand1,"$sam/$sam.txt");
my %uni=(); $count=0; my %sta=(); my %siz=(); my $number=0;
while (<hand1>)            {
	$_ =~ s/\s+$//; $count++;
	$uni{$_}++;
	$sta{substr($_,0,1)}++;
	$siz{length($_)}++;       
	}
$number=scalar keys %uni;
close hand1;
unlink "$sam/$sam.txt";
#unlink "fq/$sam.fastq";

$id=10000000;
open (hand1,">$sam/$sam\_$min\_uni.txt");
foreach my $seq (sort {$uni{$b}<=>$uni{$a}} keys %uni)  {
	$id++;
	print hand1 ">$id\_x$uni{$seq}\n$seq\n";   
	}
close hand1; %uni=(); 

open (hand1,">$sam/sum\_$sam.txt");
print hand1 "$sam\n\nTotal number of $min nt or longer reads\t$count\n";
print hand1 "Total unique species of $min nt or longer\t$number\n";
close hand1;

open (hand1,">>$sam/sum\_$sam.txt");
print hand1 "\nstart:\n";
foreach my $st (keys %sta)           {
	print hand1 $st,"\t",$sta{$st},"\n"; }
%sta=();

print hand1 "\nsize:\n";
foreach my $si (sort {$a<=>$b} keys %siz)   {
	print hand1 $si,"\t",$siz{$si},"\n";        }

close hand1; %siz=(); %sta=();              
}
