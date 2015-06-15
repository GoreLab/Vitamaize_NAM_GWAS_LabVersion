#!/usr/bin/perl

# Jason Cepela, Master Coder
# Buell Lab, Michigan State University
# 05/21/14

# This script takes a chromosome number, start position, stop position, and an output file as arguments. It reads through the HapMapv1_AGPv2, HapMapv2_AGPv2, RNA-Seq, CNV, & RDV data sets for SNPs that fall in the specified region. These data sets should be placed in a folder named "data" which should be placed in the same directory as this script.

# To run this script, type
#   perl get_me_my_SNPs_v3.pl --chr 10 --start 3828150 --stop 3828350 -o results.txt
# Please note the use of two dashes before 'chr', 'start', and 'stop' and the use of a single dash before 'o', this is important and will cause an error if not used properly.
# If you do not want the results file in your current working directory, you can use -o /path/to/file/results.txt

# Notes:
# Some lines in the CNV and RDV data files contain 0's and 1's instead of A's, C's, T's, & G's. In these instances, 0 --> A, 1 --> T.
# Some lines in the HapMap v1 & v2 files contain +'s, -'s, and 0's. In these instances, + --> A, - --> T, 0 --> W

# Version 2 updates:
# Version 1 of the script did not take into account that the RNA-seq data set excluded three samples (CML52, KI11, M162W). In version 2, results returned from the RNA-seq data set contain "NA" in these three columns.
# Version 2 of the script now outputs to the terminal any lines of input for which corresponding data is not found.

# Version 3 updates:
# This script now takes for input: chromosome number, start position, stop position, output file and finds the corresponding SNPs from the data files.

# Version 3.1 updates:
# modified by CBK to accept union files named "gwas_snps_union_chrX.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt" (where X = 1:10) used for GWAS 
#	in data_input folder in lab server /share/maize_gwas directory 
# data in this file is already major allele imputed


use warnings;
use strict;
use Getopt::Long;

my $hapmapv2_chr1 = "data_input/gwas_snps_union_chr1.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt";
my $hapmapv2_chr2 = "data_input/gwas_snps_union_chr2.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt";
my $hapmapv2_chr3 = "data_input/gwas_snps_union_chr3.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt";
my $hapmapv2_chr4 = "data_input/gwas_snps_union_chr4.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt";
my $hapmapv2_chr5 = "data_input/gwas_snps_union_chr5.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt";
my $hapmapv2_chr6 = "data_input/gwas_snps_union_chr6.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt";
my $hapmapv2_chr7 = "data_input/gwas_snps_union_chr7.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt";
my $hapmapv2_chr8 = "data_input/gwas_snps_union_chr8.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt";
my $hapmapv2_chr9 = "data_input/gwas_snps_union_chr9.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt";
my $hapmapv2_chr10 = "data_input/gwas_snps_union_chr10.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt";


# Check to ensure all data files are in place. If they are not, die.
if ( ! -e $hapmapv2_chr1 ) { die "Missing data file: $hapmapv2_chr1 \n";}
if ( ! -e $hapmapv2_chr2 ) { die "Missing data file: $hapmapv2_chr2 \n";}
if ( ! -e $hapmapv2_chr3 ) { die "Missing data file: $hapmapv2_chr3 \n";}
if ( ! -e $hapmapv2_chr4 ) { die "Missing data file: $hapmapv2_chr4 \n";}
if ( ! -e $hapmapv2_chr5 ) { die "Missing data file: $hapmapv2_chr5 \n";}
if ( ! -e $hapmapv2_chr6 ) { die "Missing data file: $hapmapv2_chr6 \n";}
if ( ! -e $hapmapv2_chr7 ) { die "Missing data file: $hapmapv2_chr7 \n";}
if ( ! -e $hapmapv2_chr8 ) { die "Missing data file: $hapmapv2_chr8 \n";}
if ( ! -e $hapmapv2_chr9 ) { die "Missing data file: $hapmapv2_chr9 \n";}
if ( ! -e $hapmapv2_chr10 ) { die "Missing data file: $hapmapv2_chr10  \n";}

my ( $output, $in_chr, $in_start, $in_stop );
my $tab = "\t";

# Get command line arguments and ensure everything is as it should be
my $usage = "\n Usage: $0 --chr <num> --start <num> --stop <num> -o <output file>";
Getopt::Long::GetOptions ( 'chr=i' => \$in_chr,
			   'start=i' => \$in_start,
			   'stop=i' => \$in_stop,
			   'o=s' => \$output );

if ( !defined($in_chr) ) { die $usage . "\n" . "Please specify a chromosome using --chr followed by a number \n"; }
if ( !defined($in_start) ) { die $usage . "\n" . "Please specify a start position using --start followed by a number \n"; }
if ( !defined($in_stop) ) { die $usage . "\n" . "Please specify a stop position using --stop followed by a number \n"; }
if ( !defined($output) ) { die $usage . "\n" . "Please specify an output file using -o followed by the name of the file \n"; }
if (   -e $output ) { print "$output already exists :( \n\n"; die; }

open ( OUT, ">", $output ) or die "Failed to open output file $output.\n$!\n";

print OUT "rs#" .$tab. "alleles" .$tab. "chrom" .$tab. "pos" .$tab. "strand" .$tab. "assembly#" .$tab. "center" .$tab. "protLSID" .$tab. "assayLSID" .$tab. "panelLSID" .$tab. "QCcode" .$tab. "B73" .$tab. "B97" .$tab. "CML103" .$tab. "CML228" .$tab. "CML247" .$tab. "CML277" .$tab. "CML322" .$tab. "CML333" .$tab. "CML52" .$tab. "CML69" .$tab. "HP301" .$tab. "IL14H" .$tab. "KI11" .$tab. "KI3" .$tab. "KY21" .$tab. "M162W" .$tab. "M37W" .$tab. "MO17" .$tab. "MO18W" .$tab. "MS71" .$tab. "NC350" .$tab. "NC358" .$tab. "OH43" .$tab. "OH7B" .$tab. "P39" .$tab. "TX303" .$tab. "TZI8" . "\n";


############################################################
print "working on HapMapv2 data set  \n";
my $hapmapv2 = "";
if ( $in_chr eq "1" ){ $hapmapv2 = $hapmapv2_chr1; } # figure out which chromsome of HapMap2 to open
if ( $in_chr eq "2" ){ $hapmapv2 = $hapmapv2_chr2; }
if ( $in_chr eq "3" ){ $hapmapv2 = $hapmapv2_chr3; }
if ( $in_chr eq "4" ){ $hapmapv2 = $hapmapv2_chr4; }
if ( $in_chr eq "5" ){ $hapmapv2 = $hapmapv2_chr5; }
if ( $in_chr eq "6" ){ $hapmapv2 = $hapmapv2_chr6; }
if ( $in_chr eq "7" ){ $hapmapv2 = $hapmapv2_chr7; }
if ( $in_chr eq "8" ){ $hapmapv2 = $hapmapv2_chr8; }
if ( $in_chr eq "9" ){ $hapmapv2 = $hapmapv2_chr9; }
if ( $in_chr eq "10" ){ $hapmapv2 = $hapmapv2_chr10; }


open ( DATA, "<", $hapmapv2 ) or die "Failed to open HapMapv2 file $hapmapv2 .\n$!\n";
<DATA>; # waster header line
while ( my $line = <DATA> )
  {
    chomp $line;                                                                             ##
    my @data = split("\t", $line);

    for my $i (11 .. $#data) # for each elemet of data                                         #
      {
        if ( $data[$i] eq "0" ) { $data[$i] = 'A'; } # if data is 0, use A                     #modified according to CNV section
		    if ( $data[$i] eq "1" ) { $data[$i] = 'T'; } # if data is 1, use T                     #modified according to CNV section
      }

    my $rs = $data[0];
    my $alleles = $data[1];
    my $chr = $data[2];
    my $pos = $data[3];
    my $strand = $data[4];
    my $assembly = $data[5];
    my $center = $data[6];
    my $protLSID = $data[7];
    my $assayLSID = $data[8];
    my $panelLSID = $data[9];
    my $QCcode = $data[10];
    my $B73 = $data[11];
    my $B97 = $data[12];
    my $CML103 = $data[13];
    my $CML228 = $data[14];
    my $CML247 = $data[15];
    my $CML277 = $data[16];
    my $CML322 = $data[17];
    my $CML333 = $data[18];
    my $CML52 = $data[19];
    my $CML69 = $data[20];
    my $HP301 = $data[21];
    my $IL14H = $data[22];
    my $KI11 = $data[23];
    my $KI3 = $data[24];
    my $KY21 = $data[25];
    my $M162W = $data[26];
    my $M37W = $data[27];
    my $MO17 = $data[28];
    my $MO18W = $data[29];
    my $MS71 = $data[30];
    my $NC350 = $data[31];
    my $NC358 = $data[32];
    my $OH43 = $data[33];
    my $OH7B = $data[34];
    my $P39 = $data[35];
    my $TX303 = $data[36];
    my $TZI8 = $data[37];
    if ( $in_start <= $pos and $pos <= $in_stop )
      {
	print OUT $rs .$tab. $alleles .$tab. $chr .$tab. $pos .$tab. $strand .$tab. $assembly .$tab. $center .$tab. $protLSID .$tab. $assayLSID .$tab. $panelLSID .$tab. $QCcode .$tab. $B73 .$tab. $B97 .$tab. $CML103 .$tab. $CML228 .$tab. $CML247 .$tab. $CML277 .$tab. $CML322 .$tab. $CML333 .$tab. $CML52 .$tab. $CML69 .$tab. $HP301 .$tab. $IL14H .$tab. $KI11 .$tab. $KI3 .$tab. $KY21 .$tab. $M162W .$tab. $M37W .$tab. $MO17 .$tab. $MO18W .$tab. $MS71 .$tab. $NC350 .$tab. $NC358 .$tab. $OH43 .$tab. $OH7B .$tab. $P39 .$tab. $TX303 .$tab. $TZI8 . "\n";
      }
  }
close DATA;

############################################################

close OUT;
print "\n\n Script has finished running. Have an awesome day! \n\n";
exit;
