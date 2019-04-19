use warnings;
use strict;

unless( @ARGV == 1 ){
	die "Wrong number of arguments.\nUsage: perl $0 <working_directory>\n";
}

my ($working_directory) = @ARGV;

my $cluster_file_1 = "$working_directory/naive_bayes_cluster1.txt";
my $cluster_file_2 = "$working_directory/naive_bayes_cluster2.txt";

my $raw_tc_cluster_pieces_file_tab = "$working_directory/../tc_sliced_vector.dat";
my $raw_tc_cluster_pieces_file_bed = "$working_directory/../tc_sliced_vector.bed";

my %clstr;
open(IN,$cluster_file_1) or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	if( exists $clstr{"1"}{$line} ){
		die "ERROR1\n";
	}

	$clstr{"1"}{$line} = 1;
# 	die "$line\n";
}
close IN;

open(IN,$cluster_file_2) or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	if( exists $clstr{"2"}{$line} ){
		die "ERROR2\n";
	}

	$clstr{"2"}{$line} = 1;
# 	die "$line\n";
}
close IN;

my %tpm;
open(IN,$raw_tc_cluster_pieces_file_tab) or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	my @temp = split(/\t/,$line);

	#!!!!!!!!Warning, the difference between replicate tpm is calculated only for rep1 and rep2
	if( exists $clstr{"1"}{$temp[0]} ){
		push @{$tpm{"1"}},abs( abs($temp[1]) - abs($temp[2]) );
# 		die abs( abs($temp[1]) - abs($temp[2]) )."\n";
	}
	elsif( exists $clstr{"2"}{$temp[0]} ){
		push @{$tpm{"2"}},abs( abs($temp[1]) - abs($temp[2]) );
# 		die abs( abs($temp[1]) - abs($temp[2]) )."\n";
	}
	else{
		die "ERROR3\n";
	}
}
close IN;

my $final_cluster;
my $cluster1_mean = mean_me(\@{$tpm{"1"}});#die "$cluster1_mean\n";
my $cluster2_mean = mean_me(\@{$tpm{"2"}});#die "$cluster2_mean\n";
if( $cluster1_mean > $cluster2_mean ){
	$final_cluster = "2";
}
elsif( $cluster2_mean > $cluster1_mean ){
	$final_cluster = "1";
}
else{
	die "ERROR4\n";
}

open(POS,">$working_directory/tc_pieces_reproducible.bed") or die "$!\n";
open(NEG,">$working_directory/tc_pieces_irreproducible.bed") or die "$!\n";
open(IN,$raw_tc_cluster_pieces_file_bed) or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	my @temp = split(/\t/,$line);

	if( exists $clstr{$final_cluster}{$temp[3]} ){
		print POS "$line\n";
	}
	else{
		print NEG "$line\n";
	}
}
close IN;
close NEG;
close POS;

############### Subroutine for calculating mean #############
sub mean_me {

	my ($array_ref) = @_;

	my $sum;
	foreach my $ele (@$array_ref){
		$sum += $ele;
	}

	return $sum/scalar(@$array_ref);
}
################################################