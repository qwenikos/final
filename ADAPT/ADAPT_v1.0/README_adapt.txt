=============================================================================
				Introduction             
=============================================================================

ADAPT is a versatile computational framework for analyzing CAGE data providing
annotation-agnostic, highly accurate and single-nucleotide resolution experimentally
derived TSSs on a genome-wide scale. It distinguishing CAGE tag-clusters that represent
transcription initiation events from clusters that are formulated due to recapping events,
byproducts of the splicing machinery as well as transcriptional and/or sequencing noise.
The sequence surrounding each CAGE tag-clusters is utilized to extract sequence and structural
DNA features which are then forwarded into a multilayered Machine Learning module.
This module is based on Support Vector (SVM) and Stochastic Gradient Boosting (SGB)
models individually trained on each DNA feature and subsequently combined with
an agent assembly strategy. Current version 1.0.

=============================================================================
				Installation         
=============================================================================

ADAPT is written in Perl (version 5 and above required) and requires the 
following Modules, all of which can be downloaded from CPAN:
 
1) Getopt::Long (it usually comes with the default perl installation)
2) File::Basename (it usually comes with the default perl installation)
3) File::Spec (it usually comes with the default perl installation)
4) PerlIO::gzip
5) Set::IntervalTree
6) Parallel::ForkManager
7) Algorithm::ExpectationMaximization;





SVM module libsvm v3.17 is required and should be in the PATH. Libsvm can
be found in http://www.csie.ntu.edu.tw/~cjlin/libsvm/.

Bedtools suite is required and should be in the PATH. Bedtools can be found
in http://code.google.com/p/bedtools/.

=============================================================================
				Using ADAPT
=============================================================================

ADAPT directory structure:

/path/to/installation/dir/ADAPT_v1.0/

	data/
		human/
			hg19/ 	#Directory that will store the hg19 results
				annotation/ #Default annotation (miRBase v20)
		mouse/
			mm9/ 	#Directory that will store the mm9 results
				annotation/ #Default annotation (miRBase v20)
			mm10/ 	#Directory that will store the mm10 results
				annotation/ #Default annotation (miRBase v20)

	example/ #Directory with data and results serving as an example

	src/ #Directory containing the source code of ADAPT
		lib/ 	#Directory containing the libraries
		model/ 	#SVM models for H3K4me3, Pol2 and DGFs
		tools/	#bam to bed conversion utility

			~~~~~~~~~~~~~~~~~~~~

Important note! You should change line 20 in ADAPT_v1.0/src/lib/algorithm_init.pm 
with the absolute ADAPT installation path in your system. I.e. line 20 should be
my $algFolder = '/home/user/ADAPT_v1.0';

ADAPT.pl is located in the installation directory of ADAPT and is the 
core script of the algorithm. By simply running ADAPT.pl without any
arguments, it produces the following output:

Usage: perl ADAPT.pl [OPTIONS] --species=<mm10/mm9/hg19> --rnaseq_directory=<RNA-seq directory> --histone_chipseq_file=<H3K4me3 ChIP-seq file> --polymerase_chipseq_file=<Pol2 ChIP-seq file> --dgf_dnase_file=<DGF DNase-seq file>

	[OPTIONS]
	--help					Print this help message.
	--file_to_scan=<String>			Name of the file containing the pre-miRNA annotation. If not set, the default annotation file will be used.
	--result_directory_id=<String>		Name of the result directory (just the name, not the full path). If left empty, a random id will be generated.
	--max_processes=<INT>			Number of processes for parallel execution of the algorithm. Default = 1 process.
	--upstream_search_region=<INT>		Upstream region to be scanned for RNA-seq enrichment. Default = 250,000 base pairs.
	--window_size=<INT>			Size of the scanning window. Default = 30 base pairs.
	--window_step=<INT>			Scanning window step in bp. Default = 5 base pairs.
	--gap_size=<INT>			Maximum required distance between discovered RNA-seq enriched islands to be consider as a single island. Default = 200.
	--window_score_threshold=<INT>		Minimum number of overlapping reads required to form an RNA-seq enriched island. Default = 5.
	--svm_score_threshold=<FLOAT>		SVM score threshold for filtering TSS predictions, 0-3. Default = 1.5.
	--in_between_prediction_cov=<FLOAT>	RNA-seq coverage threshold between TSS predictions for a given pre-miRNA. Default = 0.1.
	--strand_mode=<String>			Strand mode allowing the algorithm to consider the strandness of the RNA-seq experiment. Choose between stranded_RF/stranded_FR/unstranded. Default = unstranded.

Arguments explained:

--species 			(mandatory)	Choose between mm9,mm10 and hg19. For additional species please contact us.
--rnaseq_directory		(mandatory)	This should be the directory containing the RNA-seq reads in bed format. The
						algorithm expects to find one file per chromosome (compressed with gzip or
						not), i.e. chr13.bed.gz or chr13.bed. It is highly recommended to use the 
						bam to bed utility provided by ADAPT which is located in ADAPT_v1.0/src/tools/
--histone_chipseq_file		(mandatory)	File containing the H3K4me3 reads in bed format. It can be either compressed
						with gzip or not, i.e. H3K4me3.bed.gz or H3K4me3.bed. 
--polymerase_chipseq_file	(mandatory)	File containing the Polymerase II reads in bed format. The file can be either
						compressed with gzip or not, i.e. PolII.bed.gz or PolII.bed.
--dgf_dnase_file		(mandatory)	File containing the DGF Transcription Factor binding sites in bed format. It
						can be either compressed with gzip or not, i.e. DGF.bed.gz or DGF.bed.
--file_to_scan			(optional)	Name of the file containing the pre-miRNA annotation. If not set, the default
						annotation file will be used. The algorithm could be also used to identify
						TSS of protein coding genes and long non-coding RNAs in general. For additional
						information and instructions on this topic please contact us.
--result_directory_id		(optional)	Choose a name for the results directory (just the name, not the full path). If
						left empty, a random id will be generated. I.e. if the name "test_name" is chosen
						for mm10 species, /path/to/installation/dir/ADAPT_v1.0/data/mouse/mm10/test_name/
						will be created to hold the results.
--max_processes			(optional)	Number of processes for parallel execution of the algorithm. ADAPT has been 
						developed to work in a per chromosome fashion. Depending on the size of the 
						RNA-seq and ChIP-seq files a desktop PC with an ordinary memory size (4GB) should
						use more than 2 processes. Default = 1 process.
--upstream_search_region	(optional)	Upstream region to be scanned for RNA-seq enrichment. ADAPT applies a sliding
						window starting at the locations derived form --file_to_scan or the default miRNA
						annotation. Default = 250,000 base pairs.
--window_size			(optional)	Size of the scanning window. This option can significantly affect the results.
						It should be noted that the final prediction length will be equal to the window
						size. Default = 30 base pairs.
--window_step			(optional)	Scanning window step in bp. This should be between 1bp and the chosen window length.
						Default = 5 base pairs.
--gap_size			(optional)	Maximum required distance between two discovered RNA-seq enriched islands in order
						to be considered as a single island. This option can significantly affect the 
						results. Default = 200.
--window_score_threshold	(optional)	Minimum number of RNA-seq reads overlapping with the sliding window. This option 
						can significantly affect the results. Default = 5.
--svm_score_threshold		(optional)	SVM score threshold for filtering TSS predictions. This should be between 0-3.
						Default = 1.5.
--in_between_prediction_cov	(optional)	RNA-seq coverage threshold between TSS predictions for a given pre-miRNA. If the
						--upstream_search_region is set too large then there might be several TSS prediction
						for each pre-miRNA. This option will set the minimum RNA-seq coverage required to
						link two TSS predictions. This should be between 0-1. This option can significantly 
						affect the results. Default = 0.1.
--strand_mode			(optional)	Strand mode allowing the algorithm to consider the strandness of the RNA-seq data.
						Choose between stranded_RF (meant for Reverse-Forward paired end RNA-seq),
						stranded_FR (meant for Forward-Reverse paired end RNA-seq), unstranded (meant for
						single end RNA-seq). In the case of paired end RNA-seq data please make sure if
						it is Reverse-Forward or Forward-Reverse since this option can significantly affect
						the results. Default = unstranded.

			~~~~~~~~~~~~~~~~~~~~

ADAPT output explained:

rnaseq_enriched_islands.bed			This is the output of the first step of the algorithm which is the 
						rnaseq_enriched_region_caller.pl script. It contains all the identified RNA-seq
						enriched regions upstream of the regions included in --file_to_scan or the default
						miRNA annotation.
rnaseq_enriched_islands_scored_*.txt		This is the output of the second step of the algorithm which is the
						score_rnaseq_enriched_regions.pl script. Each RNA-seq enriched region identified
						in the previous step is scored based on H3K4me3 and Pol2 ChIP-seq and DGF TF
						binding sites as identified by DNase-seq. * refers to histone, polymerase and
						footprint.
rnaseq_enriched_islands_libsvmScored_*.txt	Same as above but this is the output of each SVM model. * refers to histone, 
						polymerase and footprint.
tss_predictions.bed				This is the output of the third step of the algorithm which is the promoter_caller.pl & 
						process_final_predictions.pl. This contains the final TSS predictions of ADAPT. Note
						that for clustered miRNAs there will be one prediction and the names will be delimited
						with ";".

			~~~~~~~~~~~~~~~~~~~~

ADAPT example run:

RNA-seq (GSM973235) and DGF (GSE40869) example data can be found in ADAPT_v1.0/example/. H3K4me3 (GSM723017) and Pol2 (GSM723019) ChIP-seq example data can be downloaded from ftp://83.212.111.252
(username: ftpuser, password: freeftp) under public/ADAPT directory and should be placed in the example directory. Output will be placed under ADAPT_v1.0/data/mouse/mm10/.

perl ADAPT.pl --file_to_scan=example/let7d_cluster.bed --result_directory_id=example_run --svm_score_threshold=1.5 --max_processes=1 --strand_mode=stranded_RF --species=mm10 --rnaseq_directory=example/ --histone_chipseq_file=example/H3K4me3.mESC.mm10.bed.gz --polymerase_chipseq_file=example/Pol2.mESC.mm10.bed.gz --dgf_dnase_file=example/DGF.mESC.mm10.bed.gz

=============================================================================
				Problems and comments 
=============================================================================

Please address any problems or comments to: 

    georgakilas@uth.gr
