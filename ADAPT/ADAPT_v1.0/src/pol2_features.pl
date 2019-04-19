use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use fur;

unless( @ARGV == 3 ){
	die "Wrong number of arguments.\nUsage: perl $0 <in_dir> <out_dir> <psem_file>\n";
}

my ($in_dir,$out_dir,$psem_file) = @ARGV;

my $fasta_file = "$in_dir/extended_tc_rep_flank100.valid.fa";
my $loc_file = "$in_dir/extended_tc_rep_flank100.valid.bed";

my %locs;
open(IN,$loc_file) or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	my @temp = split(/\t/,$line);

	@{$locs{$temp[3]}} = @temp;
}
close IN;
my $AnnotateSoftDir = "src/tools/ANNOTATE_v3.04/ANNOTATE_v3.04";
#my $AnnotateSoftDir = "/mnt/raid0/nikos/alternative_transctiption/ANNOTATE_v3.04/ANNOTATE_v3.04";
##system "/home/georgaki/bin/ANNOTATE-3.05/install_dir/bin/ANNOTATE_v3.05 -s $fasta_file --psem $psem_file --tab -o $out_dir/TRAP.tab -g 0.5";
##system "/home/nikos/bioapps/ANNOTATE_v3.05/ANNOTATE_v3.05 -s $fasta_file --psem $psem_file --tab -o $out_dir/TRAP.tab -g 0.5";
system "$AnnotateSoftDir -s $fasta_file --psem $psem_file --tab -o $out_dir/TRAP.tab -g 0.5";
open(OUT,">$out_dir/TRAP.csv") or die "$!\n";
open(OUT1,">$out_dir/TRAP.bed") or die "$!\n";
open(IN,"$out_dir/TRAP.tab") or die "$!\n";
while(my $line = <IN>){
	chomp $line;

	#command_line= /home/georgaki/bin/ANNOTATE-3.05/install_dir/bin/ANNOTATE_v3.05 -s ../extended_tc_rep.fa --psem /mnt/raid0/georgaki/projects/miRNA_genes/DATASETS/JASPAR/jaspar_pol2/cluster_buster/all_pol2.psem --tab -o TRAP.tab -g 0.5 
# 	SEQID	MTE	INR	GC-Box	CCAAT-Box	DPE	BREu	BREd	DCE_S_I	DCE_S_II	DCE_S_III	XCPE1	TATA-Box	MED1
# 	chr14_tc_31299	1.04116e-06	1.28469	2.37341e-05	0.0446017	1.22791	0.00135666	2.13642	0.236765	0.1287	0.24093	3.47147e-06	0.00238092	0.00090182

	my @temp = split(/\t/,$line);

	if( $line =~ /^#/ ){
		next;
	}
	elsif( $line =~ /^SEQID/ ){

		print OUT "label,$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$temp[8],$temp[9],$temp[10],$temp[11],$temp[12],$temp[13]\n";
	}
	else{

		if( rand(1) < 0.5 ){
			print OUT "Neg,$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$temp[8],$temp[9],$temp[10],$temp[11],$temp[12],$temp[13]\n";
		}
		else{
			print OUT "Pos,$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$temp[8],$temp[9],$temp[10],$temp[11],$temp[12],$temp[13]\n";
		}

		print OUT1 join("\t",@{$locs{$temp[0]}})."\n";
	}
}
close IN;
close OUT;
close OUT1;