use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use algorithm_init;

unless(@ARGV == 4){
	die "Wrong number of arguments.\nUsage: perl $0 <rep_file> <genome> <out_dir> <species>\n";
}

my ($rep_file,$genome,$out_dir,$species) = @ARGV;

my $window_large = 500;
my $window_large_mid = 400;

my $window_narrow = 100;
my $window_narrow_mid = 200;

system "mkdir $out_dir" unless -d "$out_dir";

my %specs = return_species_info($species);

open(LARGE,">$out_dir/extended_tc_rep_flank$window_large.bed") or die "$!\n";
open(NARROW,">$out_dir/extended_tc_rep_flank$window_narrow.bed") or die "$!\n";
open(LARGE_MID,">$out_dir/extended_tc_rep_flank$window_large_mid.bed") or die "$!\n";
open(NARROW_MID,">$out_dir/extended_tc_rep_flank$window_narrow_mid.bed") or die "$!\n";
open(IN,$rep_file) or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	my @temp = split(/\t/,$line);

	if( $temp[0] eq "chromosome" ){
		next;
	}

	my $orig_start = $temp[1];
	my $orig_stop = $temp[2];

	my $foul = 4;

	$temp[1] = $orig_start - $window_large;
	$temp[2] = $orig_stop + $window_large;
	if( $temp[1] < 0 or $temp[2] > $specs{'chr_size'}{$temp[0]} ){
		$foul--;
	}

	$temp[1] = $orig_start - $window_narrow;
	$temp[2] = $orig_stop + $window_narrow;
	if( $temp[1] < 0 or $temp[2] > $specs{'chr_size'}{$temp[0]} ){
		$foul--;
	}

	$temp[1] = $orig_start - $window_large_mid;
	$temp[2] = $orig_stop + $window_large_mid;
	if( $temp[1] < 0 or $temp[2] > $specs{'chr_size'}{$temp[0]} ){
		$foul--;
	}

	$temp[1] = $orig_start - $window_narrow_mid;
	$temp[2] = $orig_stop + $window_narrow_mid;
	if( $temp[1] < 0 or $temp[2] > $specs{'chr_size'}{$temp[0]} ){
		$foul--;
	}

	######
	if( $foul != 4 ){
		next;
	}

	$temp[1] = $orig_start - $window_large;
	$temp[2] = $orig_stop + $window_large;
	print LARGE join("\t",@temp)."\n";

	$temp[1] = $orig_start - $window_narrow;
	$temp[2] = $orig_stop + $window_narrow;
	print NARROW join("\t",@temp)."\n";

	$temp[1] = $orig_start - $window_large_mid;
	$temp[2] = $orig_stop + $window_large_mid;
	print LARGE_MID join("\t",@temp)."\n";

	$temp[1] = $orig_start - $window_narrow_mid;
	$temp[2] = $orig_stop + $window_narrow_mid;
	print NARROW_MID join("\t",@temp)."\n";
}
close IN;
close LARGE;
close NARROW;
close LARGE_MID;
close NARROW_MID;

system "bedtools getfasta -name -s -fi $genome -bed $out_dir/extended_tc_rep_flank$window_large.bed -fo $out_dir/extended_tc_rep_flank$window_large.fa";
system "bedtools getfasta -name -s -fi $genome -bed $out_dir/extended_tc_rep_flank$window_narrow.bed -fo $out_dir/extended_tc_rep_flank$window_narrow.fa";
system "bedtools getfasta -name -s -fi $genome -bed $out_dir/extended_tc_rep_flank$window_large_mid.bed -fo $out_dir/extended_tc_rep_flank$window_large_mid.fa";
system "bedtools getfasta -name -s -fi $genome -bed $out_dir/extended_tc_rep_flank$window_narrow_mid.bed -fo $out_dir/extended_tc_rep_flank$window_narrow_mid.fa";

my %invalid;

my ($header);
open(OUT,">$out_dir/extended_tc_rep_flank$window_large.tab") or die "$!\n";
open(IN,"$out_dir/extended_tc_rep_flank$window_large.fa") or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	if( $line =~ />(.+)/ ){
		$header = $1;
	}
	else{
		print OUT "$header\t".uc($line)."\n";

		if( $line =~ /N/i ){
			$invalid{$header} = 1;
		}
	}
}
close IN;
close OUT;

open(OUT,">$out_dir/extended_tc_rep_flank$window_narrow.tab") or die "$!\n";
open(IN,"$out_dir/extended_tc_rep_flank$window_narrow.fa") or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	if( $line =~ />(.+)/ ){
		$header = $1;
	}
	else{
		print OUT "$header\t".uc($line)."\n";

		if( $line =~ /N/i ){
			$invalid{$header} = 1;
		}
	}
}
close IN;
close OUT;

open(OUT,">$out_dir/extended_tc_rep_flank$window_large_mid.tab") or die "$!\n";
open(IN,"$out_dir/extended_tc_rep_flank$window_large_mid.fa") or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	if( $line =~ />(.+)/ ){
		$header = $1;
	}
	else{
		print OUT "$header\t".uc($line)."\n";

		if( $line =~ /N/i ){
			$invalid{$header} = 1;
		}
	}
}
close IN;
close OUT;

open(OUT,">$out_dir/extended_tc_rep_flank$window_narrow_mid.tab") or die "$!\n";
open(IN,"$out_dir/extended_tc_rep_flank$window_narrow_mid.fa") or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	if( $line =~ />(.+)/ ){
		$header = $1;
	}
	else{
		print OUT "$header\t".uc($line)."\n";

		if( $line =~ /N/i ){
			$invalid{$header} = 1;
		}
	}
}
close IN;
close OUT;

####### Remove 'Ns' and fix order
my @flanks = ( 'flank100', 'flank200', 'flank400', 'flank500' );
foreach my $flank (@flanks){

	open(OUT,">$out_dir/extended_tc_rep_$flank.valid.tab") or die "$!\n";
	open(OUT1,">$out_dir/extended_tc_rep_$flank.valid.fa") or die "$!\n";
	open(IN,"$out_dir/extended_tc_rep_$flank.tab") or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);

		if( exists $invalid{$temp[0]} ){
			next;
		}

		print OUT "$line\n";
		print OUT1 ">$temp[0]\n$temp[1]\n";
	}
	close IN;
	close OUT;
	close OUT1;

	open(OUT,">$out_dir/extended_tc_rep_$flank.valid.bed") or die "$!\n";
	open(IN,"$out_dir/extended_tc_rep_$flank.bed") or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);

		if( exists $invalid{$temp[3]} ){
			next;
		}

		print OUT "$line\n";
	}
	close IN;
	close OUT;
}