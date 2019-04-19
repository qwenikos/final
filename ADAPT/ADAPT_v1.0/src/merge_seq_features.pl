use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use fur;

unless( @ARGV == 2 ){
	die "Wrong number of arguments.\nUsage: perl $0 <in_dir> <out_dir>\n";
}

my ($in_dir,$out_dir) = @ARGV;

system "mkdir $out_dir" unless -d "$out_dir";

my @ft = (
	'duplex_free_energy','stacking_energy','denaturation','duplex_disrupt_energy',
	'protein_deformation','propeller_twist','z_dna','dna_bending_stiffness',
	'a_philicity','nucleosome_position','protein_dna_twist','b_dna_twist','bendability'
);

my (%combined,%locs);
foreach my $type (@ft){

	open(IN,"$in_dir/$type.scored.bed") or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);

		$combined{$temp[3]}{$type} = $temp[4];
		@{$locs{$temp[3]}} = @temp;
	}
	close IN;
}

open(OUT1,">$out_dir/combined.csv") or die "$!\n";
open(OUT2,">$out_dir/combined.bed") or die "$!\n";

print OUT1 join(",",@ft).",label\n";

foreach my $name (keys %combined){

	my @final;
	foreach my $type (@ft){

		push @final,$combined{$name}{$type};
	}

	print OUT2 join("\t",@{$locs{$name}})."\n";

	if( rand(1) < 0.5 ){

		print OUT1 join(",",@final).",Pos\n";
	}
	else{

		print OUT1 join(",",@final).",Neg\n";
	}
}
close OUT1;
close OUT2;