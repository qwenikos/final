#################################################################################################
#################################################################################################
##
##	Copyright (c) DIANA Lab, University of Thessaly, Department of Electrical and Computer Engineering, Greece
##
##	Author:	George Georgakilas
##			georgakilas@uth.gr
##
##	microTSS v2.0 - library of frequently used routines
##
##	Library containing frequently used routines by microTSS
##
#################################################################################################
#################################################################################################

use strict;
use PerlIO::gzip;
use Set::IntervalTree;
use Parallel::ForkManager;
##use Algorithm::ExpectationMaximization;

###################### Routine for finding per nucleotide CAGE occupancy ######################
sub find_CTSS {

	my ($dir,$chr) = @_;

	#die if the $file did not contain any reads for $chr
	unless(-e "$dir/$chr.bed.gz" or -e "$dir/$chr.bed"){ warn "\t\t>$dir did not contain any files for $chr\n"; }

	#create an interval tree object to load all reads of $chr
	my $tree_plus = Set::IntervalTree->new;
	my $tree_minus = Set::IntervalTree->new;
	my $check=0;

	if(-e "$dir/$chr.bed.gz"){
		open(IN,"<:gzip","$dir/$chr.bed.gz") or die "Cannot open $dir/$chr.bed.gz: $!.\n";
	}
	else{
		open(IN,"$dir/$chr.bed") or die "Cannot open $dir/$chr.bed: $!.\n";
	}
	while(my $line=<IN>){
		chomp $line;
	
		my @temp=split(/\t/, $line);

		#load $chr reads into the interval tree object
		if($temp[0] eq $chr){

			if($temp[5] eq "+"){
				$tree_plus->insert($temp[3],$temp[1],$temp[2]);
				$check=1;
			}
			elsif($temp[5] eq "-"){
				$tree_minus->insert($temp[3],$temp[1],$temp[2]);
				$check=1;
			}
		}
	}
	close IN;

	#die if the $file did not contain any reads for $chr
	unless($check == 1){ die "\t\t>>$dir/$chr.bed did not contain any reads for $chr\n"; }

	return $tree_plus,$tree_minus;
}
##########################################################################

################## Routine for converting sam to bed, one file per chromosome #################
sub convert_CAGE_sam_to_bed_per_chr {

	my ($inFile,$outDir,$cage_map_qual) = @_;

	my %valid_sam_flags = (
		'0' => '+',
		'16' => '-',
		'65' => '+',
		'67' => '+',
		'81' => '-',
		'83' => '-',
		'99' => '+',
		'115' => '-',
		'129' => '+',
		'131' => '+',
		'145' => '-',
		'147' => '-',
		'163' => '+',
		'179' => '-'
	);

	my $library_depth = 0;
	my $last_chrom = "NAN";
	my @outLine;

	open(IN,$inFile) or die "$inFile: $!\n";
	while(my $line = <IN>){
		chomp $line;

# 		VHE-1-3-6-0-33367	16	chr1	10575	3	6M1I41M	*	0	0TACGCGAGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAACTCCGCCG	*	NM:i:3	MD:Z:0C0T45	XP:Z:F$-42+"+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		my @temp = split(/\t/,$line);

		if( !exists $valid_sam_flags{$temp[1]} or $temp[4] < $cage_map_qual ){
			next;
		}

		my $tag_start = $temp[3];
		my $tag_stop = $temp[3];
		#!!!!!!!Check if chr is included in the chromosome field
		my $tag_chrom = $temp[2];
		my $tag_strand = $valid_sam_flags{$temp[1]};
		my @tag_seq = split(//,$temp[9]);
		my $tag_name = $temp[0];
		my @md_flag = split(//,$temp[12]);

		if( $tag_chrom ne $last_chrom){

			if( $last_chrom eq "NAN" ){
				$last_chrom = $tag_chrom;
			}
			else{
				open(OUT,">>$outDir/$last_chrom.bed") or die "$!\n";
				foreach my $element (@outLine){
					print OUT join("\t",@$element)."\n";
				}
				close OUT;

				$last_chrom = $tag_chrom;
				@outLine = ();
			}
		}

		my ($operators,$extensions) = split_cigar($temp[5]);

# 		if( $$operators[0] eq "I" or $$operators[0] eq "D" or $$operators[0] eq "N" ){
# 			warn "1st operator: Invalid operator $$operators[0] in CIGAR: $temp[5].\n";
# 		}
		if( $$operators[0] eq "D" or $$operators[0] eq "N" ){
			warn "1st operator: Invalid operator $$operators[0] in CIGAR: $temp[5].\n";
		}

		for(my $i = 0; $i < @$operators; $i++){

			if( ($$operators[$i] eq "I" or $$operators[$i] eq "S" or $$operators[$i] eq "H") and ($i == 0 or $i == $#{$operators}) ){
				next;
			}
			elsif( ($$operators[$i] eq "S" or $$operators[$i] eq "H") and ($i != 0 and $i != $#{$operators}) ){
				die "Invalid operator $$operators[$i] in the middle of CIGAR: $temp[5].\n";
			}
			else{
				if( $$operators[$i] eq "M" or $$operators[$i] eq "D" or $$operators[$i] eq "N" or $$operators[$i] eq "=" or $$operators[$i] eq "X" or $$operators[$i] eq "P" ){
					$tag_stop += $$extensions[$i];
				}
			}
		}

		if( $md_flag[5] == 0 and $tag_seq[0] =~ /G/i ){

			if( $tag_strand eq "+" ){
				$tag_start = $tag_start - 1 + 1;
				$tag_stop = $tag_stop - 1;
			}
			else{
				$tag_stop = $tag_stop - 2;
				$tag_start = $tag_start - 1;
			}
		}
		else{
			if( $tag_strand eq "+" ){
				$tag_start = $tag_start - 1;
				$tag_stop = $tag_stop - 1;
			}
			else{
				$tag_stop = $tag_stop - 1;
				$tag_start = $tag_start - 1;
			}
		}

		my @tmp;
		push @tmp,$tag_chrom;
		push @tmp,$tag_start;
		push @tmp,$tag_stop;
		push @tmp,$tag_name;
		push @tmp,0;
		push @tmp,$tag_strand;

		push @outLine,\@tmp;
		$library_depth++;
	}
	close IN;

	open(OUT,">>$outDir/$last_chrom.bed") or die "$!\n";
	foreach my $element (@outLine){
		print OUT join("\t",@$element)."\n";
	}
	close OUT;

	return $library_depth;
}
##########################################################################

###################### Routine for splitting CIGAR string ##########################
sub split_cigar {

	my ($string) = @_;

	my %valid_cigar_operators = (
		"M" => 1,
		"I" => 1,
		"D" => 1,
		"N" => 1,
		"S" => 1,
		"H" => 1,
		"P" => 1,
		"=" => 1,
		"X" => 1
	);

	my @temp = split(//,$string);
	my (@extensions,@operators);
	my $last_num_index = 0;
	for(my $i = 0; $i < @temp; $i++){

		unless( $temp[$i] =~ /\d/ ){

			if( !exists $valid_cigar_operators{$temp[$i]} ){
				die "Sub split_cigar: Invalid operator $temp[$i] in CIGAR: $string.\n";
			}
			else{
				push @operators,$temp[$i];
				push @extensions,join("",@temp[$last_num_index..($i-1)]);
				$last_num_index = $i + 1;
			}
		}
	}

	return \@operators,\@extensions;
}
#############################################################################

########################## Routine implementing binary search ##########################
sub binary {

	my ($array_ref,$value) = @_;

	my $min=0;
	my $max=$#{$array_ref};
	my $mid=0;

	do{
		$mid=int(($min+$max)/2);
		if($value>${${$array_ref}[$mid]}[2]){
			$min=$mid;
		}
		else{
			$max=$mid;
		}
	}until($max-$min <= 10);

	if($min-100 <=0){
		return 0;
	}
	else{
		return $min-100;
	}
}
###########################################################################

########################## Routine for loading CAGE rep file ############################
sub load_rep_file {

	my ($file) = @_;

	my %reps;
	open(IN,$file) or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);

		push @{$reps{$temp[0]}{$temp[5]}},\@temp;
	}
	close IN;

	foreach my $chr (keys %reps){
		foreach my $strand (keys %{$reps{$chr}}){
			@{$reps{$chr}{$strand}} = sort { $$a[1] <=> $$b[1] } @{$reps{$chr}{$strand}};
		}
	}

	return %reps;
}
############################################################################

########################## Routine for loading CTSS file ###############################
sub load_ctss_file {

	my ($dir,$chr) = @_;

	my %ctss;
	open(IN,"$dir/$chr.bed") or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);

		$ctss{$temp[0]}{$temp[5]}{$temp[1]} = $temp[4];
	}
	close IN;

	return %ctss;
}
############################################################################

####################### Routine for measuring dna bending stiffness #########################
sub dna_bending_stiffness {

	my ($seq) = @_;

	my %stiffness_dints = (
		'GG' => 130, 'CC' => 130,
		'GC' => 85,
		'CG' => 85,
		'AA' => 35, 'TT' => 35,
		'AT' => 20,
		'TA' => 20,
		'AC' => 60, 'GT' => 60,
		'CA' => 60, 'TG' => 60,
		'AG' => 60, 'CT' => 60,
		'GA' => 60, 'TC' => 60
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$stiffness_dints{$tmp};
	}

	return @temp;
}
############################################################################

####################### Routine for measuring dna stacking energy ##########################
sub dna_stacking_energy {

	my ($seq) = @_;

	my %stacking_energy_dints = (
		'GG' => -8.26, 'CC' => -8.26,
		'GC' => -14.59,
		'CG' => -9.69,
		'AA' => -5.37, 'TT' => -5.37,
		'AT' => -6.57,
		'TA' => -3.82,
		'AC' => -10.51, 'GT' => -10.51,
		'CA' => -6.57, 'TG' => -6.57,
		'AG' => -6.78, 'CT' => -6.78,
		'GA' => -9.81, 'TC' => -9.81
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$stacking_energy_dints{$tmp};
	}

	return @temp;
}
############################################################################

####################### Routine for measuring dna duplex free energy ########################
sub dna_duplex_free_energy {

	my ($seq) = @_;

	my %duplex_free_energy_dints = (
		'GG' => -1.77, 'CC' => -1.77,
		'GC' => -2.28,
		'CG' => -2.09,
		'AA' => -1.02, 'TT' => -1.02,
		'AT' => -0.73,
		'TA' => -0.6,
		'AC' => -1.43, 'GT' => -1.43,
		'CA' => -1.38, 'TG' => -1.38,
		'AG' => -1.16, 'CT' => -1.16,
		'GA' => -1.46, 'TC' => -1.46
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$duplex_free_energy_dints{$tmp};
	}

	return @temp;
}
############################################################################

####################### Routine for measuring protein deformation ##########################
sub protein_deformation {

	my ($seq) = @_;

	my %protein_deformation_dints = (
		'GG' => 6.1, 'CC' => 6.1,
		'GC' => 4,
		'CG' => 12.1,
		'AA' => 2.9, 'TT' => 2.9,
		'AT' => 1.6,
		'TA' => 6.3,
		'AC' => 2.3, 'GT' => 2.3,
		'CA' => 9.8, 'TG' => 9.8,
		'AG' => 2.1, 'CT' => 2.1,
		'GA' => 4.5, 'TC' => 4.5
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$protein_deformation_dints{$tmp};
	}

	return @temp;
}
############################################################################

####################### Routine for measuring protein-dna twist ###########################
sub protein_dna_twist {

	my ($seq) = @_;

	my %protein_dna_twist_dints = (
		'GG' => 32.9, 'CC' => 32.9,
		'GC' => 33.6,
		'CG' => 36.1,
		'AA' => 35.1, 'TT' => 35.1,
		'AT' => 29.3,
		'TA' => 37.8,
		'AC' => 31.5, 'GT' => 31.5,
		'CA' => 37.3, 'TG' => 37.3,
		'AG' => 31.9, 'CT' => 31.9,
		'GA' => 36.3, 'TC' => 36.3
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$protein_dna_twist_dints{$tmp};
	}

	return @temp;
}
############################################################################

######################### Routine for measuring B-dna twist #############################
sub b_dna_twist {

	my ($seq) = @_;

	my %b_dna_twist_dints = (
		'GG' => 33.4, 'CC' => 33.4,
		'GC' => 38.3,
		'CG' => 31.1,
		'AA' => 35.8, 'TT' => 35.8,
		'AT' => 33.4,
		'TA' => 40,
		'AC' => 35.8, 'GT' => 35.8,
		'CA' => 36.9, 'TG' => 36.9,
		'AG' => 30.5, 'CT' => 30.5,
		'GA' => 39.3, 'TC' => 39.3
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$b_dna_twist_dints{$tmp};
	}

	return @temp;
}
############################################################################

######################### Routine for measuring dna bendability ###########################
sub dna_bendability {

	my ($seq) = @_;

	my %dna_bendability_trints = (
		'AAA' => -0.274, 'TTT' => -0.274,
		'AAC' => -0.205, 'GTT' => -0.205,
		'AAG' => -0.081, 'CTT' => -0.081,
		'AAT' => -0.28, 'ATT' => -0.28,
		'ACA' => -0.006, 'TGT' => -0.006,
		'ACC' => -0.032, 'GGT' => -0.032,
		'ACG' => -0.033, 'CGT' => -0.033,
		'ACT' => -0.183, 'AGT' => -0.183,
		'AGA' => 0.027, 'TCT' => 0.027,
		'AGC' => 0.017, 'GCT' => 0.017,
		'AGG' => -0.057, 'CCT' => -0.057,
		'ATA' => 0.182, 'TAT' => 0.182,
		'ATC' => -0.11, 'GAT' => -0.11,
		'ATG' => 0.134, 'CAT' => 0.134,
		'CAA' => 0.015, 'TTG' => 0.015,
		'CAC' => 0.04, 'GTG' => 0.04,
		'CAG' => 0.175, 'CTG' => 0.175,
		'CCA' => -0.246, 'TGG' => -0.246,
		'CCC' => -0.012, 'GGG' => -0.012,
		'CCG' => -0.136, 'CGG' => -0.136,
		'CGA' => -0.003, 'TCG' => -0.003,
		'CGC' => -0.077, 'GCG' => -0.077,
		'CTA' => 0.09, 'TAG' => 0.09,
		'CTC' => 0.031, 'GAG' => 0.031,
		'GAA' => -0.037, 'TTC' => -0.037,
		'GAC' => -0.013, 'GTC' => -0.013,
		'GCA' => 0.076, 'TGC' => 0.076,
		'GCC' => 0.107, 'GGC' => 0.107,
		'GGA' => 0.013, 'TCC' => 0.013,
		'GTA' => 0.025, 'TAC' => 0.025,
		'TAA' => 0.068, 'TTA' => 0.068,
		'TCA' => 0.194, 'TGA' => 0.194
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 2; $i++){

		my $tmp = substr($seq,$i,3);
		push @temp,$dna_bendability_trints{$tmp};
	}

	return @temp;
}
############################################################################

###################### Routine for measuring nucleosome positioning #########################
sub nucleosome_positioning {

	my ($seq) = @_;

	my %nucleosome_positioning_trints = (
		'AAA' => -0.36, 'TTT' => -0.36,
		'AAC' => -0.06, 'GTT' => -0.06,
		'AAG' => 0.06, 'CTT' => 0.06,
		'AAT' => -0.3, 'ATT' => -0.3,
		'ACA' => 0.06, 'TGT' => 0.06,
		'ACC' => 0.08, 'GGT' => 0.08,
		'ACG' => 0.08, 'CGT' => 0.08,
		'ACT' => 0.11, 'AGT' => 0.11,
		'AGA' => -0.09, 'TCT' => -0.09,
		'AGC' => 0.25, 'GCT' => 0.25,
		'AGG' => 0.08, 'CCT' => 0.08,
		'ATA' => -0.13, 'TAT' => -0.13,
		'ATC' => 0.07, 'GAT' => 0.07,
		'ATG' => 0.18, 'CAT' => 0.18,
		'CAA' => -0.09, 'TTG' => -0.09,
		'CAC' => 0.17, 'GTG' => 0.17,
		'CAG' => -0.02, 'CTG' => -0.02,
		'CCA' => 0.08, 'TGG' => 0.08,
		'CCC' => 0.13, 'GGG' => 0.13,
		'CCG' => 0.02, 'CGG' => 0.02,
		'CGA' => 0.31, 'TCG' => 0.31,
		'CGC' => 0.25, 'GCG' => 0.25,
		'CTA' => -0.18, 'TAG' => -0.18,
		'CTC' => 0.08, 'GAG' => 0.08,
		'GAA' => -0.12, 'TTC' => -0.12,
		'GAC' => 0.08, 'GTC' => 0.08,
		'GCA' => 0.13, 'TGC' => 0.13,
		'GCC' => 0.45, 'GGC' => 0.45,
		'GGA' => -0.05, 'TCC' => -0.05,
		'GTA' => -0.06, 'TAC' => -0.06,
		'TAA' => -0.2, 'TTA' => -0.2,
		'TCA' => 0.08, 'TGA' => 0.08
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 2; $i++){

		my $tmp = substr($seq,$i,3);
		push @temp,$nucleosome_positioning_trints{$tmp};
	}

	return @temp;
}
############################################################################

######################### Routine for measuring z-dna stability############################
sub z_dna_stability {

	my ($seq) = @_;

	my %z_dna_stability_dints = (
		'GG' => 2.4, 'CC' => 2.4,
		'GC' => 4,
		'CG' => 0.7,
		'AA' => 3.9, 'TT' => 3.9,
		'AT' => 5.9,
		'TA' => 2.5,
		'AC' => 4.6, 'GT' => 4.6,
		'CA' => 1.3, 'TG' => 1.3,
		'AG' => 3.4, 'CT' => 3.4,
		'GA' => 3.4, 'TC' => 3.4
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$z_dna_stability_dints{$tmp};
	}

	return @temp;
}
############################################################################

####################### Routine for measuring duplex disrupt energy #########################
sub duplex_disrupt_energy {

	my ($seq) = @_;

	my %duplex_disrupt_energy_dints = (
		'GG' => 3.1, 'CC' => 3.1,
		'GC' => 3.1,
		'CG' => 3.6,
		'AA' => 1.9, 'TT' => 1.9,
		'AT' => 1.5,
		'TA' => 0.9,
		'AC' => 1.3, 'GT' => 1.3,
		'CA' => 1.9, 'TG' => 1.9,
		'AG' => 1.6, 'CT' => 1.6,
		'GA' => 1.6, 'TC' => 1.6
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$duplex_disrupt_energy_dints{$tmp};
	}

	return @temp;
}
############################################################################

######################### Routine for measuring propeller twist ############################
sub propeller_twist {

	my ($seq) = @_;

	my %propeller_twist_dints = (
		'GG' => -8.11, 'CC' => -8.11,
		'GC' => -11.8,
		'CG' => -10.03,
		'AA' => -18.66, 'TT' => -18.66,
		'AT' => -15.01,
		'TA' => -11.85,
		'AC' => -13.1, 'GT' => -13.1,
		'CA' => -9.45, 'TG' => -9.45,
		'AG' => -13.1, 'CT' => -14,
		'GA' => -13.48, 'TC' => -13.48
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$propeller_twist_dints{$tmp};
	}

	return @temp;
}
############################################################################

######################### Routine for measuring dna denaturation ##########################
sub dna_denaturation {

	my ($seq) = @_;

	my %dna_denaturation_dints = (
		'GG' => 118.49, 'CC' => 118.49,
		'GC' => 124.54,
		'CG' => 124.61,
		'AA' => 89.08, 'TT' => 89.08,
		'AT' => 86.72,
		'TA' => 81.85,
		'AC' => 103.18, 'GT' => 103.18,
		'CA' => 107.96, 'TG' => 107.96,
		'AG' => 99.49, 'CT' => 99.49,
		'GA' => 104.43, 'TC' => 104.43
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$dna_denaturation_dints{$tmp};
	}

	return @temp;
}
############################################################################

######################### Routine for measuring dna a-philicity ##########################
sub dna_a_philicity {

	my ($seq) = @_;

	my %dna_a_philicity_dints = (
		'GG' => 3.3, 'CC' => 3.3,
		'GC' => 3,
		'CG' => 3.9,
		'AA' => 5.7, 'TT' => 5.7,
		'AT' => 5.6,
		'TA' => 4.7,
		'AC' => 3.5, 'GT' => 3.5,
		'CA' => 3.5, 'TG' => 3.5,
		'AG' => 3.1, 'CT' => 3.1,
		'GA' => 3.9, 'TC' => 3.9
	);

	my @temp;
	for(my $i = 0; $i < length($seq) - 1; $i++){

		my $tmp = substr($seq,$i,2);
		push @temp,$dna_a_philicity_dints{$tmp};
	}

	return @temp;
}
############################################################################

########################## Routine for processing merged tc regions ########################
sub process_merged_tc_regions {

	my ($rep_working_dirs,$rep_library_depths,$working_directory,$processes) = @_;

	system "mkdir $working_directory/per_chr" unless -d "$working_directory/per_chr";

	my @rep_order;
	open(OUT,">$working_directory/tc_sliced_vector.dat") or die "$!\n";
	for(my $i = 0; $i < scalar(@$rep_working_dirs); $i++){

		$$rep_working_dirs[$i] =~ /(replicate_\d{1,2})/;
		my $rep_name = $1;
		push @rep_order,$rep_name;
	}
# 	print OUT "slice_name\t".join("\t",@rep_order)."\n";
	close OUT;

	my (%raw_spots,%uniq);
	foreach my $dirs (@$rep_working_dirs){

		warn "\t\tLoading raw tag clusters from $dirs.\n";
		open(IN,"$dirs/CAGE_tag_cluster/tc_raw/all_chr_tc.bed") or die "$!\n";
		while(my $line = <IN>){
			chomp $line;

			my @temp = split(/\t/,$line);
# 			die "$line\n";

			unless( exists $uniq{$temp[0]}{$temp[5]}{$temp[1]} ){
				push @{$raw_spots{$temp[0]}{$temp[5]}},$temp[1];
				$uniq{$temp[0]}{$temp[5]}{$temp[1]} = 1;
# 				warn "$temp[0]\t$temp[5]\t$temp[1]\n";
# 				warn "${$raw_spots{$temp[0]}{$temp[5]}}[$#{$raw_spots{$temp[0]}{$temp[5]}}]\n";
			}

			unless( exists $uniq{$temp[0]}{$temp[5]}{$temp[2]} ){
				push @{$raw_spots{$temp[0]}{$temp[5]}},$temp[2];
				$uniq{$temp[0]}{$temp[5]}{$temp[2]} = 1;
# 				warn "$temp[0]\t$temp[5]\t$temp[2]\n";
# 				warn "${$raw_spots{$temp[0]}{$temp[5]}}[$#{$raw_spots{$temp[0]}{$temp[5]}}]\n";
			}
		}
		close IN;

	}

	foreach my $chr (keys %raw_spots){
		foreach my $strand (keys %{$raw_spots{$chr}}){
# 			warn join("\n",@{$raw_spots{$chr}{$strand}})."\n";
			@{$raw_spots{$chr}{$strand}} = sort { $a <=> $b } @{$raw_spots{$chr}{$strand}};
# 			warn join("\n",@{$raw_spots{$chr}{$strand}})."\n";
		}
	}
# 	die;

	warn "\t\tLoading raw merged tag clusters from $working_directory/tc_raw.merged.sorted.bed.\n";
	my (%merged_tc);
	open(IN,"$working_directory/tc_raw.merged.sorted.bed") or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);
# 		die "$line\n$temp[0]\t$temp[4]\n";

		push @{$merged_tc{$temp[0]}{$temp[4]}},\@temp;
	}
	close IN;

	# Create a Fork Manager object to handle the child processes
	my $pm = Parallel::ForkManager->new($processes);

	foreach my $chr (keys %merged_tc){

		# Create a new child and proceed to the next one until you reach the limit of maximum number of processes
		my $pid = $pm->start and next;

		open(OUT_TAB,">$working_directory/per_chr/$chr.dat") or die "$!\n";
		open(OUT_BED,">$working_directory/per_chr/$chr.bed") or die "$!\n";

		my (%ctss,%depth);
		for(my $i = 0; $i < scalar(@$rep_working_dirs); $i++){

			$$rep_working_dirs[$i] =~ /(replicate_\d{1,2})/;
			my $rep_name = $1;

			$depth{$rep_name} = $$rep_library_depths[$i];

			open(IN,"$$rep_working_dirs[$i]/CAGE_CTSS/$chr.bed") or die "$!\n";
			while(my $line = <IN>){
				chomp $line;

				my @temp = split(/\t/,$line);
				$ctss{$rep_name}{$temp[5]}{$temp[1]} = $temp[4];
# 				die "$line\n$rep_name\t$temp[5]\t$temp[1]\t$temp[4]\n";
			}
			close IN;
		}

		foreach my $strand (keys %{$merged_tc{$chr}}){
# 			warn "$chr\t$strand\n";
			my %index = map { ${$raw_spots{$chr}{$strand}}[$_] => $_ } 0..$#{$raw_spots{$chr}{$strand}};
# 			my %index;
# 			for(my $i = 0; $i < @{$raw_spots{$chr}{$strand}}; $i++){
# 				die "${$raw_spots{$chr}{$strand}}[$i]\n";
# 			}

			for(my $i = 0; $i < @{$merged_tc{$chr}{$strand}}; $i++){

				my $start = ${${$merged_tc{$chr}{$strand}}[$i]}[1];
				my $stop = ${${$merged_tc{$chr}{$strand}}[$i]}[2];

				my $init_index = $index{$start};
				my $end_index = $index{$stop};
# 				die "$start\t$stop\n$init_index\t$end_index\n${$raw_spots{$chr}{$strand}}[$init_index]\t${$raw_spots{$chr}{$strand}}[$end_index]\n";

				for(my $j = $init_index; $j < $end_index; $j++){

					my %vector;
					for(my $k = ${$raw_spots{$chr}{$strand}}[$j]; $k < ${$raw_spots{$chr}{$strand}}[$j+1]; $k++){

						foreach my $rep_name (keys %ctss){

							if( exists $ctss{$rep_name}{$strand}{$k} ){
								$vector{$rep_name} += $ctss{$rep_name}{$strand}{$k};
							}
							else{
								$vector{$rep_name} += 0.1;
							}
						}
					}

					my @out_string;
					foreach my $rep_name (@rep_order){

						push @out_string,log_me(( ($vector{$rep_name} * 1000000) / $depth{$rep_name} ),10);
					}
					print OUT_TAB "$chr\_${$raw_spots{$chr}{$strand}}[$j]\_".${$raw_spots{$chr}{$strand}}[$j+1]."_$strand\t".join("\t",@out_string)."\n";
					print OUT_BED "$chr\t${$raw_spots{$chr}{$strand}}[$j]\t".${$raw_spots{$chr}{$strand}}[$j+1]."\t$chr\_${$raw_spots{$chr}{$strand}}[$j]\_".${$raw_spots{$chr}{$strand}}[$j+1]."_$strand\t".join("_",@out_string)."\t$strand\n";
				}
			}
		}

		close OUT_TAB;
		close OUT_BED;

		# terminates the child process
		$pm->finish;
	}

	# wait for all children and then exit
	$pm->wait_all_children;

	system "cat $working_directory/per_chr/*.dat >> $working_directory/tc_sliced_vector.dat";
	system "cat $working_directory/per_chr/*.bed > $working_directory/tc_sliced_vector.bed";
}

sub process_merged_tc_regions_v1 {

	my ($rep_working_dirs,$rep_library_depths,$working_directory,$processes) = @_;

	system "mkdir $working_directory/per_chr" unless -d "$working_directory/per_chr";

	# Load replicate order (based on user input order) and each replicate library depth
	my (@rep_order,%depth);
	for(my $i = 0; $i < scalar(@$rep_working_dirs); $i++){

		$$rep_working_dirs[$i] =~ /(replicate_\d{1,2})/;
		my $rep_name = $1;
		push @rep_order,$rep_name;

		$depth{$rep_name} = $$rep_library_depths[$i];
	}

	warn "\t\tLoading raw merged tag clusters from $working_directory/tc_raw.merged.sorted.bed.\n";
	my (%merged_tc);
	open(IN,"$working_directory/tc_raw.merged.sorted.bed") or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);
# 		die "$line\n$temp[0]\t$temp[4]\n";

		@{$merged_tc{$temp[0]}{$temp[3]}} = @temp;
	}
	close IN;

	warn "\t\tCreating replicate tpm vectors for each entry in the merged tc file.\n";
	# Create a Fork Manager object to handle the child processes
	my $pm = Parallel::ForkManager->new($processes);

	foreach my $chr (keys %merged_tc){

		# Create a new child and proceed to the next one until you reach the limit of maximum number of processes
		my $pid = $pm->start and next;

		# Load raw tag cluster tpm values per replicate, per chromosome
		my (%raw_tc_tpm);
		foreach my $dirs (@$rep_working_dirs){

			$dirs =~ /(replicate_\d{1,2})/;
			my $rep_name = $1;

			open(IN,"$dirs/CAGE_tag_cluster/tc_raw/$chr\_tc.bed") or die "$!\n";
			while(my $line = <IN>){
				chomp $line;

				my @temp = split(/\t/,$line);
# 				die "$line\n";

				$raw_tc_tpm{$rep_name}{$temp[3]} = $temp[4];
			}
			close IN;
		}

		open(OUT_TAB,">$working_directory/per_chr/$chr.dat") or die "$!\n";
		open(OUT_BED,">$working_directory/per_chr/$chr.bed") or die "$!\n";

		foreach my $tc_name (keys %{$merged_tc{$chr}}){

			my $start = ${$merged_tc{$chr}{$tc_name}}[1];
			my $stop = ${$merged_tc{$chr}{$tc_name}}[2];
			my $strand = ${$merged_tc{$chr}{$tc_name}}[4];

			my %vector;
			my @temp = split(/;/,$tc_name);

			foreach my $rep_name (@rep_order){

				my $found = 0;
				foreach my $tc_name_replicate (@temp){

					$tc_name_replicate =~ /(.+)\|(.+)/;

					if( $2 eq $rep_name ){
						$vector{$rep_name} += $raw_tc_tpm{$rep_name}{$tc_name_replicate};
						$found = 1;
					}
					else{
						next;
					}
				}

				if( $found == 0 ){
					$vector{$rep_name} = 0.1 * 1000000 / $depth{$rep_name};
				}
			}

			my @out_string;
			foreach my $rep_name (@rep_order){

				push @out_string,log_me($vector{$rep_name},10);
			}
			print OUT_TAB "$tc_name\t".join("\t",@out_string)."\n";
			print OUT_BED "$chr\t$start\t$stop\t$tc_name\t".join("_",@out_string)."\t$strand\n";

		}

		close OUT_TAB;
		close OUT_BED;

		# terminates the child process
		$pm->finish;
	}

	# wait for all children and then exit
	$pm->wait_all_children;

	system "cat $working_directory/per_chr/*.dat > $working_directory/tc_sliced_vector.dat";
	system "cat $working_directory/per_chr/*.bed > $working_directory/tc_sliced_vector.bed";
}

sub process_merged_tc_regions_v2 {

	my ($rep_working_dirs,$rep_library_depths,$working_directory,$processes) = @_;

	system "mkdir $working_directory/per_chr" unless -d "$working_directory/per_chr";

	# Load replicate order (based on user input order) and each replicate library depth
	my (@rep_order,%depth);
	for(my $i = 0; $i < scalar(@$rep_working_dirs); $i++){

		$$rep_working_dirs[$i] =~ /(replicate_\d{1,2})/;
		my $rep_name = $1;
		push @rep_order,$rep_name;

		$depth{$rep_name} = $$rep_library_depths[$i];
# 		warn "$depth{$rep_name}\n";
	}

# 	warn "\t\tLoading raw merged tag clusters from $working_directory/tc_raw.merged.sorted.bed.\n";
	warn "\t\tLoading centiled merged tag clusters from $working_directory/tc_centile.merged.sorted.bed.\n";
	my (%merged_tc);
# 	open(IN,"$working_directory/tc_raw.merged.sorted.bed") or die "$!\n";
	open(IN,"$working_directory/tc_centile.merged.sorted.bed") or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);
# 		die "$line\n$temp[0]\t$temp[4]\n";

		push @{$merged_tc{$temp[0]}{$temp[4]}},\@temp;
	}
	close IN;

	warn "\t\tCreating replicate tpm vectors for each entry in the merged tc file.\n";
	# Create a Fork Manager object to handle the child processes
	my $pm = Parallel::ForkManager->new($processes);

	foreach my $chr (keys %merged_tc){

		# Create a new child and proceed to the next one until you reach the limit of maximum number of processes
		my $pid = $pm->start and next;

		# Load raw ctss tpm values per replicate, per chromosome
		my (%ctss);
		for(my $i = 0; $i < scalar(@$rep_working_dirs); $i++){

			$$rep_working_dirs[$i] =~ /(replicate_\d{1,2})/;
			my $rep_name = $1;

			open(IN,"$$rep_working_dirs[$i]/CAGE_CTSS/$chr.bed") or die "$!\n";
			while(my $line = <IN>){
				chomp $line;

				my @temp = split(/\t/,$line);
				$ctss{$rep_name}{$chr}{$temp[5]}{$temp[1]} = $temp[4];
# 				die "$line\n$rep_name\t$temp[5]\t$temp[1]\t$temp[4]\n";
			}
			close IN;
		}

		open(OUT_TAB,">$working_directory/per_chr/$chr.dat") or die "$!\n";
		open(OUT_BED,">$working_directory/per_chr/$chr.bed") or die "$!\n";

		foreach my $strand (keys %{$merged_tc{$chr}}){

			for(my $i = 0; $i < @{$merged_tc{$chr}{$strand}}; $i++){

				my $start = ${${$merged_tc{$chr}{$strand}}[$i]}[1];
				my $stop = ${${$merged_tc{$chr}{$strand}}[$i]}[2];
				my $tc_name = ${${$merged_tc{$chr}{$strand}}[$i]}[3];

				my %vector;
				for( my $j = $start; $j < $stop; $j++){

					foreach my $rep_name (@rep_order){

						if( exists $ctss{$rep_name}{$chr}{$strand}{$j} ){
							$vector{$rep_name} += $ctss{$rep_name}{$chr}{$strand}{$j};
						}
						else{
							$vector{$rep_name} += 0;
						}
					}
				}

				my @out_string;
				foreach my $rep_name (@rep_order){

					if( $vector{$rep_name} == 0 ){
						$vector{$rep_name} = 0.1;
					}
					else{
						$vector{$rep_name} += 0.1;
					}

					push @out_string,log_me(( ($vector{$rep_name} * 1000000) / $depth{$rep_name} ),10);
				}

				print OUT_TAB "$tc_name\t".join("\t",@out_string)."\n";
				print OUT_BED "$chr\t$start\t$stop\t$tc_name\t".join("_",@out_string)."\t$strand\n";
			}
		}

		close OUT_TAB;
		close OUT_BED;

		# terminates the child process
		$pm->finish;
	}

	# wait for all children and then exit
	$pm->wait_all_children;

	system "cat $working_directory/per_chr/*.dat > $working_directory/tc_sliced_vector.dat";
	system "cat $working_directory/per_chr/*.bed > $working_directory/tc_sliced_vector.bed";
}
#############################################################################

#################### Subroutine for finding representatives in merged clusters #####################
sub find_representative_inside_merged_peaks {

	my ($rep_working_dirs,$rep_library_depths,$working_directory,$processes) = @_;

	system "mkdir $working_directory/rep_per_chr" unless -d "$working_directory/rep_per_chr";
	system "mkdir $working_directory/irrep_per_chr" unless -d "$working_directory/irrep_per_chr";

	# Load replicate order (based on user input order) and each replicate library depth
	my (@rep_order,%depth);
	for(my $i = 0; $i < scalar(@$rep_working_dirs); $i++){

		$$rep_working_dirs[$i] =~ /(replicate_\d{1,2})/;
		my $rep_name = $1;
		push @rep_order,$rep_name;

		$depth{$rep_name} = $$rep_library_depths[$i];
	}

	my %reprod;
	open(IN,"$working_directory/tc_pieces_reproducible.bed") or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);

		push @{$reprod{$temp[0]}{$temp[5]}},\@temp;
	}
	close IN;

	my %irreprod;
	open(IN,"$working_directory/tc_pieces_irreproducible.bed") or die "$!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);

		push @{$irreprod{$temp[0]}{$temp[5]}},\@temp;
	}
	close IN;

	# Create a Fork Manager object to handle the child processes
	my $pm = Parallel::ForkManager->new($processes);

	foreach my $chr (keys %reprod){

		# Create a new child and proceed to the next one until you reach the limit of maximum number of processes
		my $pid = $pm->start and next;

		# Load raw ctss tpm values per replicate, per chromosome
		my (%ctss);
		for(my $i = 0; $i < scalar(@$rep_working_dirs); $i++){

			$$rep_working_dirs[$i] =~ /(replicate_\d{1,2})/;
			my $rep_name = $1;

			open(IN,"$$rep_working_dirs[$i]/CAGE_CTSS/$chr.bed") or die "$!\n";
			while(my $line = <IN>){
				chomp $line;

				my @temp = split(/\t/,$line);
				$ctss{$rep_name}{$chr}{$temp[5]}{$temp[1]} = $temp[4];
# 				die "$line\n$rep_name\t$temp[5]\t$temp[1]\t$temp[4]\n";
			}
			close IN;
		}

		open(REP,">$working_directory/rep_per_chr/$chr.bed") or die "$!\n";

		foreach my $strand (keys %{$reprod{$chr}}){

			for(my $i = 0; $i < @{$reprod{$chr}{$strand}}; $i++){

				my $start = ${${$reprod{$chr}{$strand}}[$i]}[1];
				my $stop = ${${$reprod{$chr}{$strand}}[$i]}[2];
				my $name = ${${$reprod{$chr}{$strand}}[$i]}[3];

				my $location;
				my $max_score = 0;
				for(my $j = $start; $j < $stop; $j++){

					my $sum = 0;
					foreach my $rep_name (@rep_order){

						if( exists $ctss{$rep_name}{$chr}{$strand}{$j} ){
							$sum += ($ctss{$rep_name}{$chr}{$strand}{$j}*1000000)/$depth{$rep_name};
						}
						else{
							$sum += 0;
						}
					}

# 					if( $sum == 0 ){
# 						die "ERRORRR1\n";
# 					}

					if( $sum > $max_score ){
						$max_score = $sum;
						$location = $j;
					}
				}
				$max_score = $max_score/scalar(@rep_order);

				print REP "$chr\t$location\t".($location+1)."\t$name\t$max_score\t$strand\n";
			}
		}

		close REP;

		# terminates the child process
		$pm->finish;
	}

	# wait for all children and then exit
	$pm->wait_all_children;

	system "cat $working_directory/rep_per_chr/*.bed > $working_directory/representatives_reproducible.bed";

	# Create a Fork Manager object to handle the child processes
	$pm = Parallel::ForkManager->new($processes);

	foreach my $chr (keys %irreprod){

		# Create a new child and proceed to the next one until you reach the limit of maximum number of processes
		my $pid = $pm->start and next;

		# Load raw ctss tpm values per replicate, per chromosome
		my (%ctss);
		for(my $i = 0; $i < scalar(@$rep_working_dirs); $i++){

			$$rep_working_dirs[$i] =~ /(replicate_\d{1,2})/;
			my $rep_name = $1;

			open(IN,"$$rep_working_dirs[$i]/CAGE_CTSS/$chr.bed") or die "$!\n";
			while(my $line = <IN>){
				chomp $line;

				my @temp = split(/\t/,$line);
				$ctss{$rep_name}{$chr}{$temp[5]}{$temp[1]} = $temp[4];
# 				die "$line\n$rep_name\t$temp[5]\t$temp[1]\t$temp[4]\n";
			}
			close IN;
		}

		open(IRREP,">$working_directory/irrep_per_chr/$chr.bed") or die "$!\n";

		foreach my $strand (keys %{$irreprod{$chr}}){

			for(my $i = 0; $i < @{$irreprod{$chr}{$strand}}; $i++){

				my $start = ${${$irreprod{$chr}{$strand}}[$i]}[1];
				my $stop = ${${$irreprod{$chr}{$strand}}[$i]}[2];
				my $name = ${${$irreprod{$chr}{$strand}}[$i]}[3];

				my $location;
				my $max_score = 0;
				for(my $j = $start; $j < $stop; $j++){

					my $sum = 0;
					foreach my $rep_name (@rep_order){

						if( exists $ctss{$rep_name}{$chr}{$strand}{$j} ){
							$sum += ($ctss{$rep_name}{$chr}{$strand}{$j}*1000000)/$depth{$rep_name};
						}
						else{
							$sum += 0;
						}
					}

# 					if( $sum == 0 ){
# 						die "ERRORRR2\n";
# 					}

					if( $sum > $max_score ){
						$max_score = $sum;
						$location = $j;
					}
				}
				$max_score = $max_score/scalar(@rep_order);

				print IRREP "$chr\t$location\t".($location+1)."\t$name\t$max_score\t$strand\n";
			}
		}

		close IRREP;

		# terminates the child process
		$pm->finish;
	}

	# wait for all children and then exit
	$pm->wait_all_children;

	system "cat $working_directory/irrep_per_chr/*.bed > $working_directory/representatives_irreproducible.bed";
}

######################## Subroutine for calling the EM algorithm ############################
sub call_EM {

	my ($file,$working_directory,$rep_dirs,$em_itterations) = @_;

	my $mask = "N";
	foreach my $ele (@$rep_dirs){
		$mask .= "1";
	}
# 	die "$mask\n";

	my $clusterer = Algorithm::ExpectationMaximization -> new(
								datafile => $file,
								mask => $mask,
								K => 2,
								max_em_iterations => $em_itterations,
								seeding => 'kmeans',
								terminal_output => 0,
							);

	$clusterer -> read_data_from_file();
	$clusterer -> seed_the_clusters();
	$clusterer -> EM();
	$clusterer -> run_bayes_classifier();

	my $clusters = $clusterer -> return_disjoint_clusters();
	#!!!!!!!!!!!!!! The original module has been modified to print to a specific directory
	$clusterer -> write_naive_bayes_clusters_to_files($working_directory);
# 	$clusterer -> write_posterior_prob_clusters_above_threshold_to_files(0.2,$working_directory);

	my $new_mask;
	for(my $i = 0; $i < scalar(@$rep_dirs); $i++){
		if( $i <=3 ){
			$new_mask .= "1";
		}
		else{
			$new_mask .= "0";
		}
	}
	$clusterer -> plot_hardcopy_clusters($new_mask);
	$clusterer -> plot_hardcopy_distributions($new_mask);
# 	foreach my $index (0..@$clusters-1){
# 		print "Cluster $index (Naive Bayes): @{$clusters->[$index]}\n\n";
# 	}
}
#############################################################################

########################### Subroutine for applying log ################################
sub log_me {

	my ($number,$base) = @_;

	return log($number)/log($base);
}
############################################################################

1;