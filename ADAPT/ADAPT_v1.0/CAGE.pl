#################################################################################################
#################################################################################################
##
##	Copyright (c) DIANA Lab, University of Thessaly, Department of Electrical and Computer Engineering, Greece
##
##	Author:	George Georgakilas
##			georgakilas@uth.gr
##
##	ADAPT-CAGE v1.0
##
##	CAGE-mediated TSS identification platform
##
#################################################################################################
#################################################################################################

use warnings;
use strict;

use Getopt::Long;
use FindBin;
use File::Spec;
use lib "$FindBin::Bin/src/lib";
use algorithm_init;
use fur;

my $cage_map_qual = 10;
my $gap_size = 50;
my $result_directory_id = int(rand(10000000));
my $svm_score_threshold = 0.5;
my $max_processes = 1;
my $min_tpm = 0.5;
my $max_size = 1000;
my $em_itterations = 50;

my ($help,$species,$cage_file,$genome_file);

# Prints usage if no command line parameters are passed or there is an unknown parameter or help option is enabled
usage() if ( @ARGV < 1 or join("",@ARGV)!~/--/ or
          ! GetOptions(
		'help' => \$help,
		'species=s' => \$species,
		'result_directory_id=s' => \$result_directory_id,
		'max_processes:i' => \$max_processes,
		'gap_size:i' => \$gap_size,
		'cage_map_qual:i' => \$cage_map_qual,
		'cage_file=s' => \$cage_file,
		'genome_file=s' => \$genome_file,
		'min_tpm=f' => \$min_tpm,
		'max_size=i' => \$max_size,
		'em_itterations:i' => \$em_itterations,
		'svm_score_threshold=f' => \$svm_score_threshold
	)
or defined $help or (!defined $species or !defined $cage_file or !defined $genome_file));

# Read the directory structure for the selected species
my %species_info = return_species_info($species);
my $database_directory = $species_info{"database"};

# Resolve absolute paths of input directories
my @input_files = split(/\@/,$cage_file);
my @absolute_cage_files;
foreach my $file (@input_files){
	push @absolute_cage_files,File::Spec->rel2abs($file);
# 	warn File::Spec->rel2abs($file)."\n";
}

my $absolute_genome_file = File::Spec->rel2abs($genome_file);

# Create the working directory for the selected job
system "mkdir $database_directory/job_$result_directory_id" unless -d "$database_directory/job_$result_directory_id";
my $working_directory = "$database_directory/job_$result_directory_id";

# Display the warning message of the enabled options
warn "\n\n> ADAPT-CAGE is starting with a job id $result_directory_id.\n\n".
	"\tInput info:\n".
	"\t\tNumber of replicates - ".scalar(@input_files)."\n".
	"\t\tCAGE file(s) - ".join("\t\t\t\t\n\t\t\t\t",@absolute_cage_files)."\n".
	"\t\t$species genome file - $absolute_genome_file\n".
	"\tRequested options:\n".
	"\t\tSpecies - $species\n".
	"\t\tCage MQ - $cage_map_qual\n".
	"\t\tNumber of processes - $max_processes\n".
	"\t\tGap size - $gap_size\n".
	"\t\tMin cluster TPM - $min_tpm\n".
	"\t\tMax cluster size (bp) - $max_size\n".
	"\t\t#EM itterations - $em_itterations\n".
	"\t\tSVM score threshold - $svm_score_threshold\n".
	"\tResults will be placed in $working_directory\n\n";

# Start the algorithm
my $time;
if( scalar(@absolute_cage_files) == 1 ){

	warn "\tInitiating module for converting CAGE tags from sam to bed format.\n";
	system "mkdir $working_directory/CAGE_tags" unless -d "$working_directory/CAGE_tags";
	$time = time();
	my $library_depth = convert_CAGE_sam_to_bed_per_chr($absolute_cage_files[0],"$working_directory/CAGE_tags",$cage_map_qual);
	warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";

	warn "\tInitiating module for finding CAGE CTSSs.\n";
	system "mkdir $working_directory/CAGE_CTSS" unless -d "$working_directory/CAGE_CTSS";
	$time = time();
	system "perl $FindBin::Bin/src/ctss_caller.pl $working_directory/CAGE_tags $working_directory/CAGE_CTSS $max_processes";
	warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";

	warn "\tInitiating module for calling CAGE tag clusters.\n";
	system "mkdir $working_directory/CAGE_tag_cluster" unless -d "$working_directory/CAGE_tag_cluster";
	$time = time();
	system "perl $FindBin::Bin/src/tc_caller.pl $working_directory/CAGE_CTSS $working_directory/CAGE_tag_cluster $max_processes $gap_size $library_depth $min_tpm $max_size";
	warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";

	warn "\tGetting the sequence around CAGE tag cluster representatives.\n";
	system "mkdir $working_directory/CAGE_tag_cluster_rep_sequence" unless -d "$working_directory/CAGE_tag_cluster_rep_sequence";
	$time = time();
	system "perl $FindBin::Bin/src/sequence_window_tc_rep.pl $working_directory/CAGE_tag_cluster/tc_centile/all_chr_tc_centile_representative.bed $absolute_genome_file $working_directory/CAGE_tag_cluster_rep_sequence $species";
	warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";
}
else{

	my @rep_working_dirs;
# 	my @rep_library_depths = qw(18470997 20418169);
	my @rep_library_depths;
	my @cat_string;

	for(my $i = 0; $i < @absolute_cage_files; $i++){

		my $repl_number = $i + 1;
		warn "\tWorking with $absolute_cage_files[$i] as replicate $repl_number.\n";

		system "mkdir $working_directory/replicate_$repl_number" unless -d "$working_directory/replicate_$repl_number";
		my $new_working_directory = "$working_directory/replicate_$repl_number";
		push @rep_working_dirs,$new_working_directory;

		warn "\t\tInitiating module for converting CAGE tags from sam to bed format.\n";
		system "mkdir $new_working_directory/CAGE_tags" unless -d "$new_working_directory/CAGE_tags";
		my $time = time();
		my $library_depth = convert_CAGE_sam_to_bed_per_chr($absolute_cage_files[$i],"$new_working_directory/CAGE_tags",$cage_map_qual);
		warn "\t\t\tFinished in ".((time()-$time)/60)." minutes.\n";
		push @rep_library_depths,$library_depth;

		warn "\t\tInitiating module for finding CAGE CTSSs.\n";
		system "mkdir $new_working_directory/CAGE_CTSS" unless -d "$new_working_directory/CAGE_CTSS";
		$time = time();
		system "perl $FindBin::Bin/src/ctss_caller.pl $new_working_directory/CAGE_tags $new_working_directory/CAGE_CTSS $max_processes";
		warn "\t\t\tFinished in ".((time()-$time)/60)." minutes.\n";

		warn "\t\tInitiating module for calling CAGE tag clusters.\n";
		system "mkdir $new_working_directory/CAGE_tag_cluster" unless -d "$new_working_directory/CAGE_tag_cluster";
		$time = time();
		system "perl $FindBin::Bin/src/tc_caller.pl $new_working_directory/CAGE_CTSS $new_working_directory/CAGE_tag_cluster $max_processes $gap_size $library_depth $min_tpm $max_size";
		warn "\t\t\tFinished in ".((time()-$time)/60)." minutes.\n";

		push @cat_string,"$new_working_directory/CAGE_tag_cluster/tc_raw/all_chr_tc.bed";
		push @cat_string,"$new_working_directory/CAGE_tag_cluster/tc_centile/all_chr_tc_centile.bed";
	}

	my $final_cat = join(" ",@cat_string);
	system "mkdir $working_directory/merged_replicates" unless -d "$working_directory/merged_replicates";
# 	system "cat $final_cat | bedtools sort -i | bedtools merge -s -nms -i | bedtools sort -i > $working_directory/merged_replicates/tc_raw.merged.sorted.bed";
	system "cat $final_cat | bedtools sort -i | bedtools merge -s -nms -i | bedtools sort -i > $working_directory/merged_replicates/tc_centile.merged.sorted.bed";

	warn "\tProcessing merged tag cluster regions.\n";
	$time = time();
	process_merged_tc_regions_v2(\@rep_working_dirs,\@rep_library_depths,"$working_directory/merged_replicates",$max_processes);
	warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";

	warn "\tInitiating EM module for assessing the reproducibility of tag clusters between replicates.\n";
	system "mkdir $working_directory/merged_replicates/EM_output" unless -d "$working_directory/merged_replicates/EM_output";
	$time = time();
	call_EM("$working_directory/merged_replicates/tc_sliced_vector.dat","$working_directory/merged_replicates/EM_output",\@rep_working_dirs,$em_itterations);
	system "perl $FindBin::Bin/src/parse_em_results.pl $working_directory/merged_replicates/EM_output";
	system "mv -t $working_directory/merged_replicates/EM_output $FindBin::Bin/cluster_plot.png $FindBin::Bin/posterior_prob_plot.png";
	find_representative_inside_merged_peaks(\@rep_working_dirs,\@rep_library_depths,"$working_directory/merged_replicates/EM_output",$max_processes);
	warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";

	warn "\tGetting the sequence around CAGE tag cluster representatives.\n";
	system "mkdir $working_directory/CAGE_tag_cluster_rep_sequence" unless -d "$working_directory/CAGE_tag_cluster_rep_sequence";
	$time = time();
	system "perl $FindBin::Bin/src/sequence_window_tc_rep.pl $working_directory/merged_replicates/EM_output/representatives_reproducible.bed $absolute_genome_file $working_directory/CAGE_tag_cluster_rep_sequence $species";
	warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";
}

warn "\tInitiating sequence feature extraction and CAGE tag cluster scoring module (phase 1 - step 1).\n";
system "mkdir $working_directory/CAGE_tag_cluster_rep_scored" unless -d "$working_directory/CAGE_tag_cluster_rep_scored";
system "mkdir $working_directory/CAGE_tag_cluster_rep_scored/sequence_features" unless -d "$working_directory/CAGE_tag_cluster_rep_scored/sequence_features";
$time = time();
system "perl $FindBin::Bin/src/seq_features.pl $working_directory/CAGE_tag_cluster_rep_sequence $working_directory/CAGE_tag_cluster_rep_scored/sequence_features/individual";
warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";

warn "\tMerging sequence feature scores (phase 1 - step 2).\n";
$time = time();
system "perl $FindBin::Bin/src/merge_seq_features.pl $working_directory/CAGE_tag_cluster_rep_scored/sequence_features/individual $working_directory/CAGE_tag_cluster_rep_scored/sequence_features/merged";
system "perl $FindBin::Bin/src/caret_predict_seq_feats.pl $FindBin::Bin/src $working_directory/CAGE_tag_cluster_rep_scored/sequence_features/merged $FindBin::Bin/src/model/sequence_features/caret/model_gbm.rds";
warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";

warn "\tInitiating Pol2-associated motif affinity scanning module (phase 2).\n";
system "mkdir $working_directory/CAGE_tag_cluster_rep_scored/pol2_features" unless -d "$working_directory/CAGE_tag_cluster_rep_scored/pol2_features";
$time = time();
system "perl $FindBin::Bin/src/pol2_features.pl $working_directory/CAGE_tag_cluster_rep_sequence $working_directory/CAGE_tag_cluster_rep_scored/pol2_features $FindBin::Bin/src/model/motifs/pol2_motifs.psem";
system "perl $FindBin::Bin/src/caret_predict_pol2_feats.pl $FindBin::Bin/src $working_directory/CAGE_tag_cluster_rep_scored/pol2_features $FindBin::Bin/src/model/pol2_features/model_gbm.rds";
warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";

warn "\tInitiating final CAGE tag cluster scoring module (phase 3).\n";
$time = time();
my $rep_number = scalar(@absolute_cage_files);
warn "phase_3.1.\n";
warn $rep_number;
system "perl $FindBin::Bin/src/merge_features.pl $working_directory/CAGE_tag_cluster_rep_scored $rep_number";
warn "phase_3.2.\n";
system "perl $FindBin::Bin/src/caret_predict_combined_feats.pl $FindBin::Bin/src $working_directory/CAGE_tag_cluster_rep_scored/combined $FindBin::Bin/src/model/combined_features/model_gbm.rds";
warn "\t\tFinished in ".((time()-$time)/60)." minutes.\n";

sub usage {

	print 	"\nUnknown option: @_\n" if ( @_ );
	print	"\nUsage: perl $0 [OPTIONS] --species=<mm10/mm9/hg19/hg38> --cage_file=<CAGE-seq file(sam, \";\" separated for replicates)> --genome_file=<species_genome_file(fasta)>\n".
		"\n\t[OPTIONS]\n".
		"\t--help\tPrint this help message.\n".
		"\t--result_directory_id=<String>\tName of the result directory (just the name, not the full path). If left empty, a random id will be generated.\n".
		"\t--max_processes=<INT>\tNumber of processes for parallel execution of the algorithm. Default = 1 process.\n".
		"\t--gap_size=<INT>\tMaximum required distance between discovered CAGE CTSSs to be consider as a single peak. Default = 50.\n".
		"\t--cage_map_qual=<INT>\tCAGE reads with MQ lower than this option will be fitlered out. Default = 10 (5% mapping error probability).\n".
		"\t--min_tpm=<FLOAT>\tMinimum tpm value of accepted clusters. Default = 0.5.\n".
		"\t--max_size=<INT>\tMaximum cluster size (bp). Default = 1000.\n".
		"\t--em_itterations=<INT>\tNumber of EM itterations for finding reproducible tag clusters. Default = 50.\n".
		"\t--svm_score_threshold=<FLOAT>\tSVM score threshold for filtering TSS predictions, 0-1. Default = 0.5\n";
	exit;
}