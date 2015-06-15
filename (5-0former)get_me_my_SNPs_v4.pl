#!/usr/bin/perl

# Jason Cepela, Master Coder
# cepelaja [at] msu.edu
# Buell Lab, Plant Biology Department, Michigan State University
# 09/04/14

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

# Version 4 updates: Imputation
# This script now takes for input: chromosome number, start position, stop position, output file, and an imputation flag (-i)
# perl get_me_my_SNPs_v4.pl --chr 10 --start 3828150 --stop 3828350 -o results.txt      <-- for non-imputed results
# perl get_me_my_SNPs_v4.pl --chr 10 --start 3828150 --stop 3828350 -o results.txt -i   <-- for imputed results
#
# Including the imputation flag on the command line will have the following effects:
# All data points with no data available (N) will be imputed with the major allele for that SNP.
# A second output file will be created, with the same base name as the specified output file, with "missing.per.marker.txt" appended
# ex: -o results.txt      <-- will result in "results.txt" and "results.missing.per.marker.txt"
# This second output file will be a tab-delimitied file with four columns, as follows: rs#   chrom    pos    No.Missing
#
# Note on RNA-Seq data:
# There are three lines (CML52, KI11, M162W) which do not have data. Without the -i option, these points are filled with "N". When the -i option is present, these will be treated as any other N and will be imputed with the major allele. Because of this, all RNA-seq SNPs will have a minimum of three SNPs that were imputed, and this will be reflected in the missing.per.marker.txt file.
#
# Performance updates: since the HapMap files are ordered, once you have passed the chromosome and position you are searching for, the remainder of the file is skipped.

use warnings;
use strict;
use Getopt::Long;

# If any data files get moved or renamed, they will need to be fixed in this section.
my $hapmapv1 = "data/maizeHapMapV1_B73RefGenV2_20110309_ALL.hmp.txt";
my $hapmapv2_chr1 = "data/maizeHapMapV2_B73RefGenV2_201203028_chr1.hmp.txt";
my $hapmapv2_chr2 = "data/maizeHapMapV2_B73RefGenV2_201203028_chr2.hmp.txt";
my $hapmapv2_chr3 = "data/maizeHapMapV2_B73RefGenV2_201203028_chr3.hmp.txt";
my $hapmapv2_chr4 = "data/maizeHapMapV2_B73RefGenV2_201203028_chr4.hmp.txt";
my $hapmapv2_chr5 = "data/maizeHapMapV2_B73RefGenV2_201203028_chr5.hmp.txt";
my $hapmapv2_chr6 = "data/maizeHapMapV2_B73RefGenV2_201203028_chr6.hmp.txt";
my $hapmapv2_chr7 = "data/maizeHapMapV2_B73RefGenV2_201203028_chr7.hmp.txt";
my $hapmapv2_chr8 = "data/maizeHapMapV2_B73RefGenV2_201203028_chr8.hmp.txt";
my $hapmapv2_chr9 = "data/maizeHapMapV2_B73RefGenV2_201203028_chr9.hmp.txt";
my $hapmapv2_chr10 = "data/maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt";
my $rnaseq = "data/filtered_SNP_matrix_AGPv2_vitamaize_April_2013_Reformatted.txt";
my $cnv = "data/cnvnator_bin500_chrALL.hmp.txt";
my $rdv = "data/rdvs_win2k_chrALL_20130307.txt";

# Check to ensure all data files are in place. If they are not, die.
if ( ! -e $hapmapv1 ) { die "Missing data file: $hapmapv1 \n";}
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
if ( ! -e $rnaseq ) { die "Missing data file: $rnaseq \n";}
if ( ! -e $cnv ) { die "Missing data file: $cnv \n";}
if ( ! -e $rdv ) { die "Missing data file: $rdv \n";}

my ( $output, $in_chr, $in_start, $in_stop );
my $impute = '';
my $impute_output;
my $tab = "\t";

# Get command line arguments and ensure everything is as it should be
my $usage = "\n Usage: $0 --chr <num> --start <num> --stop <num> -o <output file> -i";
Getopt::Long::GetOptions ( 'chr=i'   => \$in_chr,
			   'start=i' => \$in_start,
			   'stop=i'  => \$in_stop,
			   'o=s'     => \$output,
			   'i'       => \$impute
			 );

if ( !defined($in_chr) ) { die $usage . "\n" . "Please specify a chromosome using --chr followed by a number \n"; }
if ( !defined($in_start) ) { die $usage . "\n" . "Please specify a start position using --start followed by a number \n"; }
if ( !defined($in_stop) ) { die $usage . "\n" . "Please specify a stop position using --stop followed by a number \n"; }
if ( !defined($output) ) { die $usage . "\n" . "Please specify an output file using -o followed by the name of the file \n"; }
if ( -e $output ) { print "$output already exists :( \n\n"; die; }

if ( $impute eq "" ) { $impute = 0; }

if ( $impute )
  {
    # create file name for imputation summary from output file name 
    $impute_output = $output;
    $impute_output =~ s/txt$//;
    $impute_output = $impute_output . "missing.per.marker.txt";
    # check to make sure imputation summary file does not already exist
    if ( -e $impute_output ) { print "$impute_output already exists :( I don't want to over-write your file.\n\n"; die; }
  }

open ( OUT, ">", $output ) or die "Failed to open output file $output.\n$!\n";
print OUT "rs#" .$tab. "alleles" .$tab. "chrom" .$tab. "pos" .$tab. "strand" .$tab. "assembly#" .$tab. "center" .$tab. "protLSID" .$tab. "assayLSID" .$tab. "panelLSID" .$tab. "QCcode" .$tab. "B73" .$tab. "B97" .$tab. "CML103" .$tab. "CML228" .$tab. "CML247" .$tab. "CML277" .$tab. "CML322" .$tab. "CML333" .$tab. "CML52" .$tab. "CML69" .$tab. "HP301" .$tab. "IL14H" .$tab. "KI11" .$tab. "KI3" .$tab. "KY21" .$tab. "M162W" .$tab. "M37W" .$tab. "MO17" .$tab. "MO18W" .$tab. "MS71" .$tab. "NC350" .$tab. "NC358" .$tab. "OH43" .$tab. "OH7B" .$tab. "P39" .$tab. "TX303" .$tab. "TZI8" . "\n";

if ( $impute ) 
  { 
    open ( IMPUTE_OUT, ">", $impute_output ) or die "Failed to open imputation summary output file $output.\n$!\n"; 
    print IMPUTE_OUT "rs#" .$tab. "chrom" .$tab. "pos" .$tab. "No.Missing" . "\n";
  }

my ( $a_count, $t_count, $c_count, $g_count, $major_allele, $num_missing, $i, $max );

############################################################
print "\n working on HapMapv1 data set \n";

open ( DATA, "<", $hapmapv1) or die "Failed to open HapMapv1 file .\n$!\n"; #open hapmap1 data file
<DATA>; # waste header line in hapmapv1 file
while ( my $line = <DATA> )  #read lines in from HapMap1 file
  {
    chomp $line;
    my ($rs, $alleles, $chr, $pos, @data) = split($tab, $line);
    if ( $in_chr == $chr and $in_start <= $pos and $pos <= $in_stop ) 
      {
	if ( $impute ) # reset allele counters
	  {
	    $a_count = 0;
	    $t_count = 0;
	    $c_count = 0;
	    $g_count = 0;
	    $num_missing = 0;
	    $max = 0;
	    $major_allele = "N";
	  }
	for $i (0 .. $#data) # for each elemet of data
	  { 
	    if ( $data[$i] eq "+" ) { $data[$i] = 'A'; } # if data is +, use A (homozygous insertion)
	    if ( $data[$i] eq "-" ) { $data[$i] = 'T'; } # if data is -, use T (homozygous deletion)
	    if ( $data[$i] eq "0" ) { $data[$i] = 'W'; } # if data is 0, use W (heterzygous)
	    if ( $impute ) # increment appropriate allele counter
	      {
		if ( $data[$i] eq "A" ) { $a_count++; }
		if ( $data[$i] eq "T" ) { $t_count++; }
		if ( $data[$i] eq "C" ) { $c_count++; }
		if ( $data[$i] eq "G" ) { $g_count++; }
	      }
	  }

	if ( $impute ) # figure major allele
	  {
	    if ( $a_count > $max )
	      {
		$max = $a_count;
		$major_allele = "A";
	      }
	    if ( $t_count > $max ) 
	      { 
		$max = $t_count;
		$major_allele = "T"; 
	      }
	    if ( $c_count > $max ) 
	      { 
		$max = $c_count;
		$major_allele = "C"; 
	      }
	    if ( $g_count > $max ) 
	      { 
		$max = $g_count;
		$major_allele = "G"; 
	      }
	  }

	print OUT $rs .$tab. $alleles .$tab. $chr .$tab. $pos; 
	
	for $i (0 .. $#data)
	  {
	    if ( $impute and $data[$i] eq "N" ) 
	      {
		$data[$i] = $major_allele;
		$num_missing++;
	      }
	    print OUT $tab . $data[$i]; 
	  }
	
	print OUT "\n";
	if ( $impute and $num_missing > 0 ) { print IMPUTE_OUT $rs .$tab. $chr .$tab. $pos .$tab. $num_missing ."\n"; }
      }
    elsif ( $chr > $in_chr ) # if you have moved past the chromosome you are looking for, skip the rest of the file
      { last; }
  }
close DATA;
      

############################################################
print " working on HapMapv2 data set  \n";
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
    chomp $line;
    my @data = split($tab, $line);
    my $pos = $data[3];
    if ( $in_start <= $pos and $pos <= $in_stop ) 
      { 
	if ( $impute ) # reset allele counters
	  {
	    $a_count = 0;
	    $t_count = 0;
	    $c_count = 0;
	    $g_count = 0;
	    $num_missing = 0;
	    $max = 0;
	    $major_allele = "N";
	  }

	my @nam_data;
	push( @nam_data, $data[11] );
	push( @nam_data, $data[12] );
	push( @nam_data, $data[42] );
	push( @nam_data, $data[47] );
	push( @nam_data, $data[48] );
	push( @nam_data, $data[49] );
	push( @nam_data, $data[51] );
	push( @nam_data, $data[53] );
	push( @nam_data, $data[61] );
	push( @nam_data, $data[63] );
	push( @nam_data, $data[69] );
	push( @nam_data, $data[70] );
	push( @nam_data, $data[71] );
	push( @nam_data, $data[72] );
	push( @nam_data, $data[73] );
	push( @nam_data, $data[74] );
	push( @nam_data, $data[75] );
	push( @nam_data, $data[76] );
	push( @nam_data, $data[77] );
	push( @nam_data, $data[78] );
	push( @nam_data, $data[79] );
	push( @nam_data, $data[80] );
	push( @nam_data, $data[81] );
	push( @nam_data, $data[82] );
	push( @nam_data, $data[84] );
	push( @nam_data, $data[105] );
	push( @nam_data, $data[106] );

	for $i (0 .. $#nam_data) # for each elemet of data
	  { 
	    if ( $nam_data[$i] eq "+" ) { $nam_data[$i] = 'A'; } # if data is +, use A (homozygous insertion)
	    if ( $nam_data[$i] eq "-" ) { $nam_data[$i] = 'T'; } # if data is -, use T (homozygous deletion)
	    if ( $nam_data[$i] eq "0" ) { $nam_data[$i] = 'W'; } # if data is 0, use W (heterzygous)
	    if ( $impute ) # increment appropriate allele counter
	      {
		if ( $nam_data[$i] eq "A" ) { $a_count++; }
		if ( $nam_data[$i] eq "T" ) { $t_count++; }
		if ( $nam_data[$i] eq "C" ) { $c_count++; }
		if ( $nam_data[$i] eq "G" ) { $g_count++; }
	      }
	  }

	if ( $impute ) # figure major allele, change Ns to Major 
	  {
	    if ( $a_count > $max )
	      {
		$max = $a_count;
		$major_allele = "A";
	      }
	    if ( $t_count > $max ) 
	      { 
		$max = $t_count;
		$major_allele = "T"; 
	      }
	    if ( $c_count > $max ) 
	      { 
		$max = $c_count;
		$major_allele = "C"; 
	      }
	    if ( $g_count > $max ) 
	      { 
		$max = $g_count;
		$major_allele = "G"; 
	      }

	    for $i (0 .. $#nam_data) # for each elemet of data
	      { 
		if ( $nam_data[$i] eq "N" ) 
		  {  
		    $nam_data[$i] = $major_allele;
		    $num_missing++;
		  }
	      }
	  } 
		
	my $rs = $data[0];
	my $alleles = $data[1];
	my $chr = $data[2];	
	my $strand = $data[4];
	my $assembly = $data[5];
	my $center = $data[6];
	my $protLSID = $data[7];
	my $assayLSID = $data[8];
	my $panelLSID = $data[9];
	my $QCcode = $data[10];
	my $B73 = $nam_data[0];
	my $B97 = $nam_data[1];
	my $CML103 = $nam_data[2];
	my $CML228 = $nam_data[3];
	my $CML247 = $nam_data[4];
	my $CML277 = $nam_data[5];
	my $CML322 = $nam_data[6];
	my $CML333 = $nam_data[7];
	my $CML52 = $nam_data[8];
	my $CML69 = $nam_data[9];
	my $HP301 = $nam_data[10];
	my $IL14H = $nam_data[11];
	my $KI11 = $nam_data[12];
	my $KI3 = $nam_data[13];
	my $KY21 = $nam_data[14];
	my $M162W = $nam_data[15];
	my $M37W = $nam_data[16];
	my $MO17 = $nam_data[17];
	my $MO18W = $nam_data[18];
	my $MS71 = $nam_data[19];
	my $NC350 = $nam_data[20];
	my $NC358 = $nam_data[21];
	my $OH43 = $nam_data[22];
	my $OH7B = $nam_data[23];
	my $P39 = $nam_data[24];
	my $TX303 = $nam_data[25];
	my $TZI8 = $nam_data[26];
	
	print OUT $rs .$tab. $alleles .$tab. $chr .$tab. $pos .$tab. $strand .$tab. $assembly .$tab. $center .$tab. $protLSID .$tab. $assayLSID .$tab. $panelLSID .$tab. $QCcode .$tab. $B73 .$tab. $B97 .$tab. $CML103 .$tab. $CML228 .$tab. $CML247 .$tab. $CML277 .$tab. $CML322 .$tab. $CML333 .$tab. $CML52 .$tab. $CML69 .$tab. $HP301 .$tab. $IL14H .$tab. $KI11 .$tab. $KI3 .$tab. $KY21 .$tab. $M162W .$tab. $M37W .$tab. $MO17 .$tab. $MO18W .$tab. $MS71 .$tab. $NC350 .$tab. $NC358 .$tab. $OH43 .$tab. $OH7B .$tab. $P39 .$tab. $TX303 .$tab. $TZI8 . "\n";

	if ( $impute and $num_missing > 0 ) { print IMPUTE_OUT $rs .$tab. $chr .$tab. $pos .$tab. $num_missing ."\n"; }

      }
    elsif ( $in_stop < $pos ) # Since the files are ordered, if the position of the SNP is > the indicated range, stop going through the file
      { last; }
  }
close DATA;

############################################################
print " working on RDV data set \n";

open ( DATA, "<", $rdv ) or die "Failed to open RDV file.\n$!\n";
while ( my $line = <DATA> ) 
  {
    chomp $line;
    my ($rs, $alleles, $chr, $pos, @data) = split($tab, $line);
    if ( $in_chr eq $chr and $in_start <= $pos and $pos <= $in_stop ) 
      {
	
	if ( $impute ) # reset allele counters
	  {
	    $a_count = 0;
	    $t_count = 0;
	    $c_count = 0;
	    $g_count = 0;
	    $num_missing = 0;
	    $max = 0;
	    $major_allele = "N";
	  }
	
	for $i (0 .. $#data) # for each elemet of data
	  { 
	    if ( $data[$i] eq "0" ) { $data[$i] = 'A'; } # if data is 0, use A
	    if ( $data[$i] eq "1" ) { $data[$i] = 'T'; } # if data is 1, use T
	    if ( $impute ) # increment appropriate allele counter
	      {
		if ( $data[$i] eq "A" ) { $a_count++; }
		if ( $data[$i] eq "T" ) { $t_count++; }
		if ( $data[$i] eq "C" ) { $c_count++; }
		if ( $data[$i] eq "G" ) { $g_count++; }
	      }
	  }

	if ( $impute ) # figure major allele
	  {
	    if ( $a_count > $max )
	      {
		$max = $a_count;
		$major_allele = "A";
	      }
	    if ( $t_count > $max ) 
	      { 
		$max = $t_count;
		$major_allele = "T"; 
	      }
	    if ( $c_count > $max ) 
	      { 
		$max = $c_count;
		$major_allele = "C"; 
	      }
	    if ( $g_count > $max ) 
	      { 
		$max = $g_count;
		$major_allele = "G"; 
	      }
	  } 

	print OUT $rs .$tab. $alleles .$tab. $chr .$tab. $pos; 
	
	for $i (0 .. $#data)
	  { 
	    if ( $impute and $data[$i] eq "N" ) 
	      {  
		$data[$i] = $major_allele;
		$num_missing++;
	      }
	    print OUT $tab . $data[$i]; 
	  }
	
	print OUT "\n";
	if ( $impute and $num_missing > 0 ) { print IMPUTE_OUT $rs .$tab. $chr .$tab. $pos .$tab. $num_missing ."\n"; }
      }
  }

close DATA;
   
############################################################
print " working on CNV data set \n";
    	
open ( DATA, "<", $cnv ) or die "Failed to open CNV file.\n$!\n";
while ( my $line = <DATA> ) 
  {
    chomp $line;
    my ($rs, $alleles, $chr, $pos, @data) = split($tab, $line);
    if ( $in_chr eq $chr and $in_start <= $pos and $pos <= $in_stop ) 
      {
	if ( $impute ) # reset allele counters
	  {
	    $a_count = 0;
	    $t_count = 0;
	    $c_count = 0;
	    $g_count = 0;
	    $num_missing = 0;
	    $max = 0;
	    $major_allele = "N";
	  }
	
	for $i (0 .. $#data) # for each elemet of data
	  { 
	    if ( $data[$i] eq "0" ) { $data[$i] = 'A'; } # if data is 0, use A
	    if ( $data[$i] eq "1" ) { $data[$i] = 'T'; } # if data is 1, use T
	    if ( $impute ) # increment appropriate allele counter
	      {
		if ( $data[$i] eq "A" ) { $a_count++; }
		if ( $data[$i] eq "T" ) { $t_count++; }
		if ( $data[$i] eq "C" ) { $c_count++; }
		if ( $data[$i] eq "G" ) { $g_count++; }
	      }
	  }
	
	if ( $impute ) # figure major allele
	  {
	    if ( $a_count > $max )
	      {
		$max = $a_count;
		$major_allele = "A";
	      }
	    if ( $t_count > $max ) 
	      { 
		$max = $t_count;
		$major_allele = "T"; 
	      }
	    if ( $c_count > $max ) 
	      { 
		$max = $c_count;
		$major_allele = "C"; 
	      }
	    if ( $g_count > $max ) 
	      { 
		$max = $g_count;
		$major_allele = "G"; 
	      }
	  } 

	print OUT $rs .$tab. $alleles .$tab. $chr .$tab. $pos; 
	
	for $i (0 .. $#data)
	  { 
	    if ( $impute and $data[$i] eq "N" ) 
	      {  
		$data[$i] = $major_allele;
		$num_missing++;
	      }
	    print OUT $tab . $data[$i]; 
	  }
	
	print OUT "\n";
	if ( $impute and $num_missing > 0 ) { print IMPUTE_OUT $rs .$tab. $chr .$tab. $pos .$tab. $num_missing ."\n"; }
      }
  }

close DATA;
    
############################################################
print " working on RNA-seq data set \n";

open ( DATA, "<", $rnaseq ) or die "Failed to open RNA-Seq file.\n$!\n";
while ( my $line = <DATA> ) 
  {
    chomp $line;
    my @data = split($tab, $line); 
    my $chr = $data[2];
    my $pos = $data[3];
    if ( $in_chr eq $chr and $in_start <= $pos and $pos <= $in_stop )       
      { 

	if ( $impute ) 
	  {
	    # reset allele counters
	    $a_count = 0;
	    $t_count = 0;
	    $c_count = 0;
	    $g_count = 0;
	    $num_missing = 0;
	    $max = 0;
	    $major_allele = "N";

	    for $i (11 .. $#data) # count alleles
	      { 
		if ( $data[$i] eq "A" ) { $a_count++; }
		if ( $data[$i] eq "T" ) { $t_count++; }
		if ( $data[$i] eq "C" ) { $c_count++; }
		if ( $data[$i] eq "G" ) { $g_count++; }
	      }
	    
	    # figure major allele
	    if ( $a_count > $max )
	      {
		$max = $a_count;
		$major_allele = "A";
	      }
	    if ( $t_count > $max ) 
	      { 
		$max = $t_count;
		$major_allele = "T"; 
	      }
	    if ( $c_count > $max ) 
	      { 
		$max = $c_count;
		$major_allele = "C"; 
	      }
	    if ( $g_count > $max ) 
	      { 
		$max = $g_count;
		$major_allele = "G"; 
	      }

	    # replace Ns with major allele
	    for $i (11 .. $#data)
	      { 
		if ( $data[$i] eq "N" ) 
		  {  
		    $data[$i] = $major_allele;
		    $num_missing++;
		  }
	      }
	  }

	my $rs = $data[0];
	my $alleles = $data[1];
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
	my $CML52 = "N";
	my $CML69 = $data[19];
	my $HP301 = $data[20];
	my $IL14H = $data[21];
	my $KI11 = "N";
	my $KI3 = $data[22];
	my $KY21 = $data[23];
	my $M162W = "N";
	my $M37W = $data[24];
	my $MO17 = $data[26];
	my $MO18W = $data[27];
	my $MS71 = $data[25];
	my $NC350 = $data[28];
	my $NC358 = $data[29];
	my $OH43 = $data[30];
	my $OH7B = $data[31];
	my $P39 = $data[32];
	my $TX303 = $data[33];
	my $TZI8 = $data[34];

	if ( $impute ) # impute 3 missing data points, reflect this missing data in the summary output.
	  {
	    $CML52 = $major_allele;
	    $num_missing++;
	    $KI11 = $major_allele;
	    $num_missing++;
	    $M162W = $major_allele;
	    $num_missing++;
	    print IMPUTE_OUT $rs .$tab. $chr .$tab. $pos .$tab. $num_missing ."\n";
	  }

	print OUT $rs .$tab. $alleles .$tab. $chr .$tab. $pos .$tab. $strand .$tab. $assembly .$tab. $center .$tab. $protLSID .$tab. $assayLSID .$tab. $panelLSID .$tab. $QCcode .$tab. $B73 .$tab. $B97 .$tab. $CML103 .$tab. $CML228 .$tab. $CML247 .$tab. $CML277 .$tab. $CML322 .$tab. $CML333 .$tab. $CML52 .$tab. $CML69 .$tab. $HP301 .$tab. $IL14H .$tab. $KI11 .$tab. $KI3 .$tab. $KY21 .$tab. $M162W .$tab. $M37W .$tab. $MO17 .$tab. $MO18W .$tab. $MS71 .$tab. $NC350 .$tab. $NC358 .$tab. $OH43 .$tab. $OH7B .$tab. $P39 .$tab. $TX303 .$tab. $TZI8 . "\n";

      }
  }

close DATA;	
 
############################################################

close OUT;
if ( $impute ) { close IMPUTE_OUT; }
print "\n Script has finished running. Have an awesome day! \n\n";
exit;
