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

my %seq_feats = (
	'duplex_free_energy' => \&dna_duplex_free_energy,
	'stacking_energy' => \&dna_stacking_energy,
	'denaturation' => \&dna_denaturation,
	'duplex_disrupt_energy' => \&duplex_disrupt_energy,
	'protein_deformation' => \&protein_deformation,
	'propeller_twist' => \&propeller_twist,
	'z_dna' => \&z_dna_stability,
	'dna_bending_stiffness' => \&dna_bending_stiffness,
	'a_philicity' => \&dna_a_philicity,
	'nucleosome_position' => \&nucleosome_positioning,
	'protein_dna_twist' => \&protein_dna_twist,
	'b_dna_twist' => \&b_dna_twist,
	'bendability' => \&dna_bendability
);

my %optimal_windows = (
	'duplex_free_energy' => {'window' => 'flank500','size' => 'size_1'},
	'stacking_energy' => {'window' => 'flank500','size' => 'size_10'},
	'denaturation' => {'window' => 'flank500','size' => 'size_1'},
	'duplex_disrupt_energy' => {'window' => 'flank400','size' => 'size_1'},
	'protein_deformation' => {'window' => 'flank500','size' => 'size_10'},
	'propeller_twist' => {'window' => 'flank100','size' => 'size_1'},
	'z_dna' => {'window' => 'flank400','size' => 'size_10'},
	'dna_bending_stiffness' => {'window' => 'flank500','size' => 'size_10'},
	'a_philicity' => {'window' => 'flank500','size' => 'size_10'},
	'nucleosome_position' => {'window' => 'flank400','size' => 'size_1'},
	'protein_dna_twist' => {'window' => 'flank200','size' => 'size_1'},
	'b_dna_twist' => {'window' => 'flank500','size' => 'size_10'},
	'bendability' => {'window' => 'flank500','size' => 'size_1'},
);

my $rep_file = "$in_dir/extended_tc_rep_flank100.valid.bed";

foreach my $type (keys %optimal_windows){

	my $window = $optimal_windows{$type}{'window'};
	my $size = $optimal_windows{$type}{'size'};

	warn "\t\t\tWorking with $type.\n";

	open(OUT,">$out_dir/$type.tab") or die "$!\n";
	open(IN,"$in_dir/extended_tc_rep_$window.valid.tab") or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

# 		tc_name	sequence

		my @temp = split(/\t/,$line);

		$temp[1] = uc($temp[1]);

		if( $temp[1] =~ /N/ ){
			die "$temp[1]\n";
		}

		my @intermediate = $seq_feats{$type}->($temp[1]);

		if( $size eq "size_10" ){

			my @final;
			for(my $k = 0; $k < scalar(@intermediate) - 10 + 1; $k++){
				my $av;
				for(my $n = $k; $n < $k + 10; $n++){
					$av += $intermediate[$n];
				}
				push @final,$av/10;
			}

			for(my $i = 0; $i < @final; $i++){
				$final[$i] = ($i + 1).":".$final[$i];
			}

			print OUT "1 ".join(" ",@final)."\n";
		}
		else{

			for(my $i = 0; $i < @intermediate; $i++){
				$intermediate[$i] = ($i + 1).":".$intermediate[$i];
			}

			print OUT "1 ".join(" ",@intermediate)."\n";
		}
	}
	close IN;
	close OUT;
    my $SvmScaleSoftDir = "src/tools/svm/svm-scale";
	my $SvmPredictSoftDir = "src/tools/svm/svm-predict";
	system "$SvmScaleSoftDir -l -1 -u 1 -r $FindBin::Bin/model/sequence_features/libsvm/$type.scapa $out_dir/$type.tab > $out_dir/$type.scaled";
	
	my $pred =  `$SvmPredictSoftDir -b 1 $out_dir/$type.scaled $FindBin::Bin/model/sequence_features/libsvm/$type.model $out_dir/$type.scored`;
	system "perl $FindBin::Bin/translate_libsvm_output.pl $out_dir/$type.scored $rep_file > $out_dir/$type.scored.bed";

}