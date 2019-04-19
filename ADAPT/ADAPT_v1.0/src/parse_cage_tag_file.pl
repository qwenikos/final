use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use algorithm_init;
use fur;

unless( @ARGV == 2 ){
	die "Wrong number of arguments.\nUsage: perl $0 <input_sam> <out_dir>\n";
}

my ($input_sam,$out_dir) = @ARGV;

my $time = time();

warn "\tReading sam file.\n";
convert_CAGE_sam_to_bed_per_chr($input_sam,$out_dir);
warn "\tFinished in ".(time()-$time)." seconds.\n";