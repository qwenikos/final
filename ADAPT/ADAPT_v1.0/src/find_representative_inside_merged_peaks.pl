use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use fur;

unless( @ARGV == 4 ){
	die "Wrong number of arguments.\nUsage: perl $0 <rep_working_dirs> <rep_library_depths> <working_directory> <processes>\n";
}

my ($rep_working_dirs,$rep_library_depths,$working_directory,$processes) = @ARGV;

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
open(IN,"$working_directory/EM_output/tc_pieces_reproducible.bed") or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	my @temp = split(/\t/,$line);

	push @{$reprod{$temp[0]}{$temp[5]},\@temp;
}
close IN;

my %irreprod;
open(IN,"$working_directory/EM_output/tc_pieces_irreproducible.bed") or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	my @temp = split(/\t/,$line);

	push @{$irreprod{$temp[0]}{$temp[5]},\@temp;
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
# 			die "$line\n$rep_name\t$temp[5]\t$temp[1]\t$temp[4]\n";
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
						$sum += $ctss{$rep_name}{$chr}{$strand}{$j};
					}
					else{
						$sum += 0;
					}
				}

				if( $sum == 0 ){
					die "ERRORRR\n";
				}

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
# 			die "$line\n$rep_name\t$temp[5]\t$temp[1]\t$temp[4]\n";
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
						$sum += $ctss{$rep_name}{$chr}{$strand}{$j};
					}
					else{
						$sum += 0;
					}
				}

				if( $sum == 0 ){
					die "ERRORRR2\n";
				}

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