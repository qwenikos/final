use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use fur;

unless( @ARGV == 3 ){
	die "Wrong number of arguments.\nUsage: perl $0 <cage_tag_dir> <out_dir> <processes>\n";
}

my ($cage_tag_dir,$out_dir,$processes) = @ARGV;

my @files = glob("$cage_tag_dir/*.bed");

# Create a Fork Manager object to handle the child processes
my $pm = Parallel::ForkManager->new($processes);

foreach my $file (@files){

	# Create a new child and proceed to the next one until you reach the limit of maximum number of processes
	my $pid = $pm->start and next;

	$file =~ /$cage_tag_dir\/(.+)\.bed/;
	my $chr = $1;
	warn "\t\t\tWorking on $chr.\n";

	my $ctss_counter = 0;
	my %ctss;
	open(IN,$file) or die "$file: $!\n";
	while(my $line = <IN>){
		chomp $line;

		my @temp = split(/\t/,$line);

		if( $temp[5] eq "+" ){
			$ctss{$temp[5]}{$temp[1]}++;
		}
		else{
			$ctss{$temp[5]}{($temp[2] - 1)}++;
		}
	}
	close IN;

	open(OUT,">$out_dir/$chr.bed") or die "$!\n";
	foreach my $strand (keys %ctss){
		foreach my $coordinate ( sort { $a <=> $b } keys(%{$ctss{$strand}}) ){

			$ctss_counter++;

			print OUT "$chr\t$coordinate\t".($coordinate + 1)."\t$chr\_ctss_$ctss_counter\t$ctss{$strand}{$coordinate}\t$strand\n";
		}
	}
	close OUT;

	# terminates the child process
	$pm->finish;
}

# wait for all children and then exit
$pm->wait_all_children;