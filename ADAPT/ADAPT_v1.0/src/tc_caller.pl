use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use fur;

unless( @ARGV == 7 ){
	die "Wrong number of arguments.\nUsage: perl $0 <cage_ctss_dir> <out_dir> <processes> <gap_size> <library_depth> <min_tpm> <max_size>\n";
}

my ($cage_ctss_dir,$out_dir,$processes,$gap_size,$library_depth,$min_tpm,$max_size) = @ARGV;

system "mkdir $out_dir/tc_raw" unless -d "$out_dir/tc_raw";
system "mkdir $out_dir/tc_centile" unless -d "$out_dir/tc_centile";

my @files = glob("$cage_ctss_dir/*.bed");

# Create a Fork Manager object to handle the child processes
my $pm = Parallel::ForkManager->new($processes);

foreach my $file (@files){

	# Create a new child and proceed to the next one until you reach the limit of maximum number of processes
	my $pid = $pm->start and next;

	$file =~ /$cage_ctss_dir\/(.+)\.bed/;
	my $chr = $1;
	warn "\t\t\tWorking on $chr.\n";

	my (%ctss,%ctss_indexed);
	open(IN,$file) or die "$file: $!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);

		push @{$ctss{$temp[5]}},\@temp;
		$ctss_indexed{$temp[5]}{$temp[1]} = $temp[4];
	}
	close IN;

	foreach my $strand (keys %ctss){
		@{$ctss{$strand}} = sort { $$a[1] <=> $$b[1] } @{$ctss{$strand}};
	}

	my $tc_counter = 0;
	my %tc_raw;
	my $last_cluster_index = 0;
	foreach my $strand (keys %ctss){
		for(my $i = 0;  $i < @{$ctss{$strand}}; $i++){

			my $ctss_start = ${${$ctss{$strand}}[$i]}[1];
			my $ctss_stop = ${${$ctss{$strand}}[$i]}[2];
			my $ctss_score = ${${$ctss{$strand}}[$i]}[4];

			if( $i == 0 ){

				my @tmp;
				$tc_counter++;

				push @tmp,$chr;
				push @tmp,$ctss_start;
				push @tmp,$ctss_stop;
				push @tmp,"$chr\_tc_$tc_counter";
				push @tmp,$ctss_score;
				push @tmp,$strand;

				push @{$tc_raw{$strand}},\@tmp;

				$last_cluster_index = 0;
			}
			else{
				if( $ctss_start - ${${$tc_raw{$strand}}[$last_cluster_index]}[2] <= $gap_size ){
					${${$tc_raw{$strand}}[$last_cluster_index]}[2] = $ctss_stop;
					${${$tc_raw{$strand}}[$last_cluster_index]}[4] += $ctss_score;
				}
				else{
					my @tmp;
					$tc_counter++;

					push @tmp,$chr;
					push @tmp,$ctss_start;
					push @tmp,$ctss_stop;
					push @tmp,"$chr\_tc_$tc_counter";
					push @tmp,$ctss_score;
					push @tmp,$strand;

					push @{$tc_raw{$strand}},\@tmp;

					$last_cluster_index++;
				}
			}
		}
	}

	open(TC_BED,">$out_dir/tc_raw/$chr\_tc.bed") or die "$out_dir/tc_raw/$chr\_tc.bed: $!\n";
	open(TC_TAB,">$out_dir/tc_raw/$chr\_tc.tab") or die "$out_dir/tc_raw/$chr\_tc.tab: $!\n";
	my %tc_centile;
	foreach my $strand (keys %tc_raw){
		for(my $i = 0; $i < @{$tc_raw{$strand}}; $i++){

			my $tc_start = ${${$tc_raw{$strand}}[$i]}[1];
			my $tc_stop = ${${$tc_raw{$strand}}[$i]}[2];
			${${$tc_raw{$strand}}[$i]}[4] = ${${$tc_raw{$strand}}[$i]}[4] * 1000000 / $library_depth;
			my $tc_tpm = ${${$tc_raw{$strand}}[$i]}[4];

			if( $tc_stop - $tc_start > $max_size ){
				next;
			}

			#!!!!!!!!!!!Add command line argument for tpm threshold
			if( $tc_tpm < $min_tpm ){
				next;
			}

			print TC_BED join("\t",@{${$tc_raw{$strand}}[$i]})."\n";

			my $rep_start = 0;
			my $rep_score = 0;

			for(my $j = $tc_start; $j < $tc_stop; $j++){

				unless( exists $ctss_indexed{$strand}{$j} ){
					next;
				}

				if( $ctss_indexed{$strand}{$j} > $rep_score ){
					$rep_start = $j;
					$rep_score = $ctss_indexed{$strand}{$j};
				}
			}

			print TC_TAB join("\t",@{${$tc_raw{$strand}}[$i]})."\t$rep_start\t".($rep_start + 1)."\t".($rep_score * 1000000 / $library_depth)."\n";

			my ($left_start,$right_start);
			my $left_cumulative_tpm = 0;
			my $right_cumulative_tpm = 0;

			for(my $j = $tc_start; $j < $tc_stop; $j++){

				unless( exists $ctss_indexed{$strand}{$j} ){
					next;
				}

				$left_cumulative_tpm += $ctss_indexed{$strand}{$j} * 1000000 / $library_depth;
				if( $left_cumulative_tpm / $tc_tpm < 0.1 ){
					next;
				}
				else{
					$left_start = $j;
					last;
				}
			}

			for(my $j = $tc_stop - 1; $j >= $tc_start; $j--){

				unless( exists $ctss_indexed{$strand}{$j} ){
					next;
				}

				$right_cumulative_tpm += $ctss_indexed{$strand}{$j} * 1000000 / $library_depth;
				if( $right_cumulative_tpm / $tc_tpm < 0.1 ){
					next;
				}
				else{
					$right_start = $j;
					last;
				}
			}

			if( $right_start < $left_start ){
				die "ERROR1: ".join("\t",@{${$tc_raw{$strand}}[$i]})."\n$right_start\t$left_start\n";
			}
			else{

				my $tc_centile_tpm;
				for(my $j = $left_start; $j <= $right_start; $j++){

					unless( exists $ctss_indexed{$strand}{$j} ){
						next;
					}

					$tc_centile_tpm += $ctss_indexed{$strand}{$j} * 1000000 / $library_depth;
				}

				my @tmp;
				push @tmp,$chr;
				push @tmp,$left_start;
				push @tmp,$right_start + 1;
				push @tmp,${${$tc_raw{$strand}}[$i]}[3];
				push @tmp,$tc_centile_tpm;
				push @tmp,$strand;

				push @{$tc_centile{$strand}},\@tmp;
			}
		}
	}
	close TC_BED;
	close TC_TAB;

	open(CENTILE_BED,">$out_dir/tc_centile/$chr\_tc_centile.bed") or die "$out_dir/tc_centile/$chr\_tc_centile.bed: $!\n";
	open(CENTILE_TAB,">$out_dir/tc_centile/$chr\_tc_centile.tab") or die "$out_dir/tc_centile/$chr\_tc_centile.tab: $!\n";
	open(REP,">$out_dir/tc_centile/$chr\_tc_centile_representative.bed") or die "$out_dir/tc_centile/$chr\_tc_centile_representative.bed: $!\n";
	foreach my $strand (keys %tc_centile){
		for(my $i = 0; $i < @{$tc_centile{$strand}}; $i++){

			my $tc_start = ${${$tc_centile{$strand}}[$i]}[1];
			my $tc_stop = ${${$tc_centile{$strand}}[$i]}[2];
			my $tc_tpm = ${${$tc_centile{$strand}}[$i]}[4];

# 			if( $tc_tpm < 1 ){
# 				next;
# 			}

			print CENTILE_BED join("\t",@{${$tc_centile{$strand}}[$i]})."\n";

			my $rep_start = 0;
			my $rep_score = 0;

			for(my $j = $tc_start; $j < $tc_stop; $j++){

				unless( exists $ctss_indexed{$strand}{$j} ){
					next;
				}

				if( $ctss_indexed{$strand}{$j} > $rep_score ){
					$rep_start = $j;
					$rep_score = $ctss_indexed{$strand}{$j};
				}
			}

			print CENTILE_TAB join("\t",@{${$tc_centile{$strand}}[$i]})."\t$rep_start\t".($rep_start + 1)."\t".($rep_score * 1000000 / $library_depth)."\n";
			print REP "$chr\t$rep_start\t".($rep_start + 1)."\t${${$tc_centile{$strand}}[$i]}[3]\t".($rep_score * 1000000 / $library_depth)."\t$strand\n";
		}
	}
	close CENTILE_BED;
	close CENTILE_TAB;
	close REP;

	# terminates the child process
	$pm->finish;
}

# wait for all children and then exit
$pm->wait_all_children;

system "echo \"chrom\ttc_start\ttc_stop\ttc_name\ttc_tag_count_tpm\tstrand\ttc_representative_start\ttc_representative_stop\ttc_representative_score\" > $out_dir/tc_raw/all_chr_tc.tab";
system "cat $out_dir/tc_raw/chr*_tc.tab >> $out_dir/tc_raw/all_chr_tc.tab";
system "cat $out_dir/tc_raw/chr*_tc.bed > $out_dir/tc_raw/all_chr_tc.bed";
system "echo  \"chrom\ttc_percentile_start\ttc_percentile_stop\ttc_percentile_name\ttc_percentile_tag_count_tpm\tstrand\ttc_percentile_representative_start\ttc_percentile_representative_stop\ttc_percentile_representative_score\" > $out_dir/tc_centile/all_chr_tc_centile.tab";
system "cat $out_dir/tc_centile/chr*_tc_centile.tab >> $out_dir/tc_centile/all_chr_tc_centile.tab";
system "cat $out_dir/tc_centile/*_tc_centile_representative.bed > $out_dir/tc_centile/all_chr_tc_centile_representative.bed";
system "cat $out_dir/tc_centile/*_tc_centile.bed > $out_dir/tc_centile/all_chr_tc_centile.bed";