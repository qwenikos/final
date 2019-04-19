use warnings;
use strict;

unless( @ARGV == 3 ){
	die "Wrong number of arguments.\nUsage: perl $0 <source_dir> <in_dir> <model_file>\n";
}

my ($source_dir,$in_dir,$model_file) = @ARGV;

my $annotation_file = "$in_dir/representatives.bed";
my $seq_feats_file = "$in_dir/representatives.csv";

my $comm = `Rscript $source_dir/run_caret_on_combined_feats.R $model_file $seq_feats_file > $in_dir/caret_output.txt`;

my $counter = 0;
my %locs;
open(IN,$annotation_file) or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	my @temp = split(/\t/,$line);

	$counter++;
	@{$locs{$counter}} = @temp;
}
close IN;

$counter = 0;
open(OUT,">$in_dir/representatives.scored.bed") or die "$!\n";
open(IN,"$in_dir/caret_output.txt") or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	if( $line =~ /Pos/i ){
		next;
	}

	$counter++;
	my @temp = split(/\s+/,$line);

	${$locs{$counter}}[4] = $temp[2];
	print OUT join("\t",@{$locs{$counter}})."\n";
}
close IN;
close OUT;