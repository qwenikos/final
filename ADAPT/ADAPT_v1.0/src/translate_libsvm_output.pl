use warnings;
use strict;

unless(@ARGV==2){
	die "Wrong number of arguments.\nUsage: perl $0 <svmOut> <preSvm>\n";
}

my ($svmOut,$preSvm) = @ARGV;

my %svmPre;
my $counter=0;
open(IN,$preSvm) or die "$!\n";
while(my $line=<IN>){
	chomp $line;
	
	$counter++;
	my @temp=split(/\t/,$line);
	@{$svmPre{$counter}}=@temp;
}
close IN;

$counter=0;
open(IN,$svmOut) or die "$!\n";
while(my $line=<IN>){
	chomp $line;
	
	if($line=~/label/){
		next;
	}
	
	$counter++;
	
	my @temp=split(/\s/,$line);

	${$svmPre{$counter}}[4] = $temp[1];
	print join("\t",@{$svmPre{$counter}})."\n";
}
close IN;