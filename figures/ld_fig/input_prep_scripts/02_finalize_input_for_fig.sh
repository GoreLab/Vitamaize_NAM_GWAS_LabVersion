#!/bin/sh


ld_file_default='region_snps-ld_peak_snp.tsv';
expr_file_default='ExprEffEst_Correl-allQTL.tsv';
expr_signif_file_default='Signif_ExprEffEst_Correl-allQTL.tsv';
rmip_file_default='allTraits_RMIP.tsv';
gff_file_source_default='ftp://ftp.gramene.org/pub/gramene/maizesequence.org/release-5b/working-set/ZmB73_5a_WGS.gff.gz';
bedtools_bin_default="$(which bedtools)"

if [ ${#@} -le 0 ]
then
    echo "Usage: bash $0 input_directory output_directory [ld_file] [expr_file] [expr_signif_file] [rmip_file] [gff3_url]
REQUIRED:
  input_directory: directory containing combined expression data, ld data, and gwas data
  output_directory: directory where to save data figures output
OPTIONAL:
  ld_file: filename for ld input file
               default = '${ld_file_default}'
  expr_file: filename for expression correlation file
               default = '${expr_file_default}'
  expr_signif_file: filename for expression correlation pvalues file
               default = '${expr_signif_file_default}'
  rmip_file: filename for RMIP file
               default = '${rmip_file_default}'
  gff3_url: URL for the gff3 file containing genome annotation for the reference genome
               default = '${gff_file_source_default}'
  bedtools_bin: path to bedtools binary
               default = '${bedtools_bin_default}'
NOTE: Specifying any of the optional parameters requires specifying all of them in the given order
";
    exit;
fi



input_data_dir=$1
fig_data_dir=$2
ld_file=$3;
expr_file=$4;
expr_signif_file=$5;
rmip_file=$6;
gff_url=$7;
bedtools_bin=$8;

in_data_d="$(readlink -fnq "$input_data_dir")"
if [ -z "$in_data_d" ]
then
    echo "ERROR: Could not find input data directory $input_data_dir" >&2;
    exit 1
fi

fig_data_d="$(readlink -fnq "$fig_data_dir")"
if [ -z "$fig_data_d" ]
then
    echo "ERROR: Could not find output data directory $fig_data_dir" >&2;
    exit 1
fi

if [ -z "$ld_file" ]
then
    ld_file="$in_data_d/$ld_file_default"
fi
ld_f="$(readlink -fnq "$ld_file")"
if [ -z "$ld_f" ]
then
    echo "ERROR: Could not find input file ${ld_file}" >&2;
    exit 1
fi


if [ -z "$expr_file" ]
then
    expr_file="$in_data_d/$expr_file_default"
fi
expr_f="$(readlink -fnq "$expr_file")"
if [ -z "$expr_f" ]
then
    echo "ERROR: Could not find input file ${expr_file}" >&2;
    exit 1
fi

if [ -z "$expr_signif_file" ]
then
    expr_signif_file="$in_data_d/$expr_signif_file_default"
fi
expr_signif_f="$(readlink -fnq "$expr_signif_file")"
if [ -z "$expr_signif_f" ]
then
    echo "ERROR: Could not find input file ${expr_signif_file}" >&2;
    exit 1
fi


if [ -z "$rmip_file" ]
then
    rmip_file="$in_data_d/$rmip_file_default"
fi
rmip_f="$(readlink -fnq "$rmip_file")"
if [ -z "$rmip_f" ]
then
    echo "ERROR: Could not find input file ${rmip_file}" >&2;
    exit 1
fi

if [ -z "$gff_url" ]
then
    gff_url="$gff_file_source_default"
fi
gff_file="${in_data_d}/$(basename "$gff_url")"

if [ ! -f "$gff_file" ]
then  
    echo "Retrieving $gff_url"
    wget -q -O "$gff_file" "$gff_url"
fi

gff_f="$(readlink -fnq "$gff_file")"
if [ -z "$gff_f" ]
then
    echo "ERROR: Could not find local copy of ${gff_file} or retrieve using url ${gff_url}" >&2;
    exit 1
fi

if [ -z "$bedtools_bin" ]
then
    bedtools_bin="$bedtools_bin_default"
fi
bedtools_bin_f="$(readlink -fnq "$bedtools_bin")"
if [ -z "$bedtools_bin_f" ]
then 
    echo "ERROR: Could not find binary for bedtools: $bedtools_bin" >&2;
    exit 1
fi

echo "
######################
# PARAMETERS
#=====================
# INPUT DIR: $in_data_d
# OUTPUT DIR: $fig_data_d
# LD FILE: $ld_f
# EXPR FILE: $expr_f
# EXPR_SIG FILE: $expr_signif_f
# RMIP FILE: $rmip_f
# GFF_URL: $gff_url
# GFF FILE: $gff_f
# BEDTOOLS BIN: $bedtools_bin_f
######################
"

# Create list of QTLs to use from the QTL expr signif file

#==================
# NOTE: for custom QTL filtering, replace this output (qtl_filter_list.tsv) with the desired subset
#==================

qtl_filter_file="${fig_data_d}/tmp_file-qtl_filter_list.tsv"

echo "# Creating list of QTLs to keep and saving to $qtl_filter_file"

if [ -f "$qtl_filter_file" ]
then
    echo "$qtl_filter_file already exists, skipping"
else
    echo -e "QTL_No\tQTL\tGRMZM_ID" > "$qtl_filter_file"
    tail -n+2 "$expr_signif_f" | \
	cut -f1,2,4 | \
	sort -u | \
	sort -k1n,1 -k2V,2 >> "$qtl_filter_file"
fi



##########################################
# Create overall filter file with QTL, candidate gene, and peak marker, one entry per QTL
# Steps:
#  1. Find coordinates for all candidate genes
#  2. Find all the peak SNPs within or nearest each candidate gene
#  3. If multiple peak SNPs are within or equidistant to the gene, select the one with the highest overall RMIP, regardless of trait


# extract coordinates for the candidate genes
gene_coord_file="${fig_data_d}/fig_data-candidate_gene_coord.bed"

echo "# Extracting gene coordinates to $gene_coord_file"

if [ -f "$gene_coord_file" ]
then
    echo "$gene_coord_file already exists, skipping"
else
    echo -e "#chr\tstart\tstop\tGRMZM_ID\tQTL\tstrand" > "${gene_coord_file}"

    # NOTE: this grep approach works because the list is short and "clean"
    gid_list=$(cut -f3 "$qtl_filter_file" | paste -s -d'|')
    gzip -kcd "$gff_f" | \
	grep -P "Name=(${gid_list});" | \
	cut -f1,4,5,7,9- | \
	sed -re 's/\t[^\t]*Name=([^\;]+)(\;|$)[^\t]*$/\t\1/' | \
	sort -k5b,5 -t$'\t' | \
	join -1 5 -2 3 -t$'\t' - <(tail -n+2 "$qtl_filter_file" | sort -k3b,3 -t$'\t') -o 1.1,1.2,1.3,1.5,2.2,1.4 | \
	sort -k1,1 -k2,2n >> "${gene_coord_file}"
fi


# Find SNP peak marker proximal to candidate gene

all_nearest_peaks_file="${fig_data_d}/tmp_file-all_nearest_peaks.bed"
echo "# Saving all nearest peak markers for each candidate gene to $all_nearest_peaks_file"

if [ -f "$all_nearest_peaks_file" ]
then
    echo "$all_nearest_peaks_file already exists, skipping"
else

    tail -n+2 "$ld_f" | \
	cut -f5-7 | \
	sort -u | \
	awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $2, $2, $3}' | \
	sort -k1,1 -k2,2n | \
	"$bedtools_bin_f" closest -d -a "$gene_coord_file" -b stdin | \
	awk 'BEGIN{FS="\t";OFS="\t"}{print $7, $8, $9, $4":"$10, $11}' > "$all_nearest_peaks_file"
fi        


# Use RMIP as tiebreaker for equidistant SNPs
maxRMIP_nearest_peaks_file="${fig_data_d}/tmp_file-maxRMIP_nearest_peaks.bed"
echo "# Saving all selected nearest peak markers for each candidate gene to $maxRMIP_nearest_peaks_file"

if [ -f "$maxRMIP_nearest_peaks_file" ]
then
    echo "$maxRMIP_nearest_peaks_file already exists, skipping"
else

    sed -re 's/\:/\t/' -e 's/^([^\t]+)\t([^\t]+)/\1.\2\t\1\t\2/' "$all_nearest_peaks_file" | \
	sort -k1b,1 | \
	cut -f1,5 | \
	join -j1 -t $'\t' - <(tail -n+2 "$rmip_f" | sed -re 's/^([^\t]+)\t([^\t]+)/\1.\2\t\1\t\2/' | cut -f1,5,6 | sort -k1b,1) | \
	sort -k2,2 -k1,1 | \
	"$bedtools_bin_f" groupby -g 2,1 -c 3 -o max | \
	sort -k1,1 -k3nr,3 | \
	"$bedtools_bin_f" groupby -g 1 -c 2 -o first | \
	cut -f2 | \
	sort -k1b,1 | \
	join -j1 -t$'\t' - <(sed -r -e 's/^([^\t]+)\t([^\t]+)/\1.\2\t\1\t\2/' "$all_nearest_peaks_file" | sort -k1b,1) | \
	cut -f2- > "$maxRMIP_nearest_peaks_file"


fi

#=======================
# NOTE: Change this to the appropriate version of the 1peak/Gene file as needed
peaks_file="$maxRMIP_nearest_peaks_file"
#=======================

# Collate selected keys into one file
filter_file="${fig_data_d}/fig_data-combined_filter_list.tsv"

echo "# Creating list of QTLs, Genes, and Peaks to keep and saving to $filter_file"

if [ -f "$filter_file" ]
then
    echo "$filter_file already exists, skipping"
else

    (echo -e 'QTL_No\tQTL\tGRMZM_ID\tPeakChr\tPeakPos\tPeakAllele' && \
	awk 'BEGIN{FS="\t|:"; OFS="\t"}{print $4,$1,$2,$5}' "$peaks_file" | \
	sort -k1b,1 | \
	join -1 1 -2 3 -t$'\t' - <(tail -n+2 "$qtl_filter_file" | sort -k3b,3) -o 2.1,2.2,2.3,1.2,1.3,1.4 -a 2 ) > "$filter_file"
fi


#####################################################


# Filter gene expression file and add gene names
gene_expr_corr_file="${fig_data_d}/fig_data-gene_expr_correl.tsv"

echo "# Filtering gene expression correlation data, adding gene id, and saving to $gene_expr_corr_file"

if [ -f "$gene_expr_corr_file" ]
then 
    echo "$gene_expr_corr_file already exists, skipping"
else
    (head -n1 "$expr_f" | \
	sed -re 's/\tTrait/\tGRMZM_ID\tTrait/' && \
	tail -n+2 "$filter_file" | \
	awk 'BEGIN{FS="\t";OFS="\t"}{print $1"."$2, $1, $2, $3}' | \
	sort -k1b,1 | \
	join -t$'\t' -j 1 - <(tail -n+2 "$expr_f" | sed -re 's/^([^\t]+)\t([^\t]+)/\1.\2/' | sort -k1b,1) | \
	cut -f2-) > "$gene_expr_corr_file" 
fi


# Filter gene expression correlation significance file
gene_expr_corr_sig_file="${fig_data_d}/fig_data-gene_expr_correl_sig.tsv"

echo "# Filtering gene expression correlation significance data and saving to $gene_expr_corr_sig_file"

if [ -f "$gene_expr_corr_sig_file" ]
then 
    echo "$gene_expr_corr_sig_file already exists, skipping"
else
    (head -n1 "$expr_signif_f" && \
	tail -n+2 "$filter_file" | \
	awk 'BEGIN{FS="\t";OFS="\t"}{print $1"."$2}' | \
	sort -k1b,1 | \
	join -t$'\t' -j 1 - <(tail -n+2 "$expr_signif_f" | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"."$2, $0}' | sort -k1b,1) | \
	cut -f2-) > "$gene_expr_corr_sig_file" 
fi


# filter ld file based on selected peak markers
ld_filtered_file="${fig_data_d}/fig_data-ld_from_peaks.tsv"
echo "# Filtering ld data for kept peak snps and saving to $ld_filtered_file"

if [ -f "$ld_filtered_file" ] 
then
    echo "$ld_filtered_file already exists, skipping"
else

    ( head -n1 "$ld_f" && \
	cut -f4,5 "$filter_file" | \
	sed -re 's/\t/\./' | \
	sort -k1b,1 | \
	join -j1 -t$'\t' - <(tail -n+2 "$ld_f" | awk 'BEGIN{FS="\t";OFS="\t"}{print $5"."$6, $0}' | sort -k1b,1) | \
	cut -f2- ) > "$ld_filtered_file"

fi



# filter rmip file based on QTL_No ID
rmip_filtered_file="${fig_data_d}/fig_data-gwas_rmip.tsv"
echo "# Filtering gwas hits for kept QTLs and saving to $rmip_filtered_file"

if [ -f "$rmip_filtered_file" ]
then
    echo "$rmip_filtered_file already exists, skipping"
else

    qtlNo_list=$(tail -n+2 "$filter_file" | cut -f1 | sort -u | paste -s -d'|'); 
    (head -n1 "$rmip_f" && \
	tail -n+2 "$rmip_f" | \
	grep -P "\t(${qtlNo_list})\s*$" ) > "$rmip_filtered_file"
fi
