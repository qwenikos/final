use warnings;
use strict;

unless( @ARGV == 2 ){
	die "Wrong number of arguments.\nUsage: perl $0 <in_dir> <rep_number>\n";
}

my ($in_dir,$rep_number) = @ARGV;

system "mkdir $in_dir/combined" unless -d "$in_dir/combined";

my %original;
if( $rep_number == 1 ){
	open(IN,"$in_dir/../CAGE_tag_cluster/tc_centile/all_chr_tc_centile_representative.bed") or die "$!\n";
	while( my $line = <IN> ){
		chomp $line;

		my @temp = split(/\t/,$line);

		$temp[4] = log10($temp[4]);

		@{$original{$temp[3]}} = @temp;
	}
	close IN;
}
else{
	open(IN,"$in_dir/../merged_replicates/EM_output/representatives_reproducible.bed") or die "$!\n";
	while( my $line = <IN> ){
		chomp $line;

		my @temp = split(/\t/,$line);

		$temp[4] = log10($temp[4]);

		@{$original{$temp[3]}} = @temp;
	}
	close IN;
}

my $seq_feats_file = "$in_dir/sequence_features/merged/combined.scored.bed";
my $pol2_feats_file = "$in_dir/pol2_features/TRAP.scored.bed";

my %score;
open(IN,$seq_feats_file) or die "$seq_feats_file: $!\n";
while(my $line = <IN>){
	chomp $line;

	my @temp = split(/\t/,$line);

	$score{$temp[3]}{"SEQ_FEATS"} = $temp[4];
}
close IN;

open(IN,$pol2_feats_file) or die "$pol2_feats_file: $!\n";
while(my $line = <IN>){
	chomp $line;

	my @temp = split(/\t/,$line);

	$score{$temp[3]}{"POL2_FEATS"} = $temp[4];
}
close IN;

open(OUT1,">$in_dir/combined/representatives.bed") or die "$in_dir/combined/representatives.bed: $!\n";
open(OUT2,">$in_dir/combined/representatives.csv") or die "$in_dir/combined/representatives.csv: $!\n";

print OUT2 "label,EXPR,SEQ_FEATS,POL2_FEATS\n";

foreach my $name (keys %score){

	if( rand(1) < 0.5 ){

		print OUT2 "Pos,${$original{$name}}[4],".$score{$name}{"SEQ_FEATS"}.",".$score{$name}{"POL2_FEATS"}."\n";
	}
	else{

		print OUT2 "Neg,${$original{$name}}[4],".$score{$name}{"SEQ_FEATS"}.",".$score{$name}{"POL2_FEATS"}."\n";
	}

	print OUT1 join("\t",@{$original{$name}})."\n";
}
close OUT1;
close OUT2;

sub log10 {

	my ($in_num) = @_;

	return log($in_num) / log(10);
}