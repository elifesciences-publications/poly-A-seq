use strict;
use warnings;

# Feb 14, 2017
# Yunkun Dang


# Bin Tian's method, 4 random nucleotides before T stretch
# trim the PASS reads based on Bin Tian's Nature Method 2013.
# set the minimum length of T as 1.
## the sequence is reversed. 


# set minimum number of T
my $min=1;
# set minimum length
my $size = 20;


my @sample = ("SRR2225336", "SRR2225337", "SRR2225338", "SRR2225339", "SRR2225340", "SRR2225341", "SRR2225342", "SRR2225343", "SRR2225344", "SRR2225345", "SRR2225346", "SRR2225347", "SRR2225348", "SRR2225349", "SRR2225350", "SRR2225351");

foreach my $id (@sample)   {
	print "processing $id now\n\n";
	trim3($id);  }



####trim 3'
sub trim3                                                      {
my ($sam)=@_; my $count=0; my $id=10000000; my $idfq; my $on=0; my $si=0;
my $cca=0;my %cca=(); my $t = 0;

mkdir $sam;
mkdir 'fq';
open (hand1,"fq/$sam.fastq") or die $!;
open (hand2,">$sam/$sam\_A_$min.txt");
open (hand3,">fq/$sam\_A_$min.fq");

while (<hand1>)                 {
	$_ =~ s/\s+$//;
	$_ = substr($_, 0, 50);
	$count++;
	my $mod=$count % 4;
	my $len = length($_);

	if ($mod == 1)       {
		$_ =~ s/\s+.+$//g;
		$idfq=$_; next  }

	if ($mod == 2)             {
		$on=0; $si=0;
	if (/^([ATCG]{4}T{$min,}.+)$/)   {
		$_=$1;
		$_ =~ s/^[ATCG]{4}T+//;   }
	
	next if length($_) < $size;
	next if /N/;
	
	if ($len - length($_) >= 4 + $min)   {
		$on=1; $id++;
		$si=length($_);
		$t = $len-$si-4;
		print hand2 $_,"\n";
		print hand3 "$idfq\_A=$t\n";
		my $seq = reverse $_;
		$seq =~ tr/ATGC/TACG/; 
		print hand3 $seq,"\n";
		print hand3 "+$id\n";}     }
	elsif ($mod ==0 && $on==1) {
		my $qua=reverse(substr($_,$t+4,$si));
		print hand3 $qua,"\n";     }   }
close hand1; close hand2; close hand3; 


#### uni
open (hand1,"$sam/$sam\_A_$min.txt");
my %uni=(); $count=0; my %sta=(); my %siz=(); my $number=0;
while (<hand1>)            {
$_ =~ s/\s+$//; $count++;
$uni{$_}++;
$sta{substr($_,0,1)}++;
$siz{length($_)}++;        }
$number=scalar keys %uni;
close hand1;
unlink "$sam/$sam\_A_$min.txt";
#unlink "fq/$sam.fastq";

$id=10000000;
open (hand1,">$sam/$sam\_$min\_uni.txt");
foreach my $seq (sort {$uni{$b}<=>$uni{$a}} keys %uni)  {
$id++;
print hand1 ">$id\_x$uni{$seq}\n$seq\n";                }
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

close hand1; %siz=(); %sta=();                                   }
