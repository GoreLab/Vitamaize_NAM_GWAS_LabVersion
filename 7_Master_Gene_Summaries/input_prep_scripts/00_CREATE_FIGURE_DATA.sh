#!/bin/sh

# NOTE: component scripts are expected in the same directory as current script
combine_script="$(dirname $0)/01_prepare_input.sh"
filter_script="$(dirname $0)/02_finalize_input_for_fig.sh"

combined_dir_name_default='combined'

ld_file_tag_default='.target_ld.'
expr_file_tag_default='Expr.EffEst_Correl_QTL'
expr_signif_file_tag_default='Signif_of_Expr.EffEst.Correls'
rmip_file_tag_default='_RMIPgt4_with_QTLnumber_forFig2'

ld_file_default='region_snps-ld_peak_snp.tsv';
expr_file_default='ExprEffEst_Correl-allQTL.tsv';
expr_signif_file_default='Signif_ExprEffEst_Correl-allQTL.tsv';
rmip_file_default='allTraits_RMIP.tsv';
gff_file_source_default='ftp://ftp.gramene.org/pub/gramene/maizesequence.org/release-5b/working-set/ZmB73_5a_WGS.gff.gz';
bedtools_bin_default="$(which bedtools)"


if [ ${#@} -le 0 ]
then
    echo "Usage: bash $0 original_data_directory figure_data_directory [combined_data_directory]
REQUIRED:
  original_data_directory: directory containing expression data, ld data, and gwas data subdirectories
  figure_data_directory: directory where to save transformed data that will be used as input to plotting
OPTIONAL:
  combined_data_directory: an intermediate files directory where processed and combined original data will be saved before filtering for figure (default = subdirectory '${combined_dir_name_default}' in the original data directory)
";
    exit;
fi

in_dir=$1; #directory containing expression data, ld data, and gwas data subdirectories
out_dir=$2; #location of data to use as input for the figure
comb_dir=$3; 

in_d="$(readlink -fnq "$in_dir")"
if [ -z $in_d ]
then
    echo "ERROR: Could not find input directory $in_dir" >&2;
    exit 1;
fi

out_d="$(readlink -fnq "$out_dir")"
if [ -z $out_d ]
then
    echo "ERROR: Could not find output directory $out_dir" >&2;
    exit 1;
fi

if [ -z $comb_dir ]
then
    comb_dir="${in_d}/${combined_dir_name_default}"
fi
[ -d "$comb_dir" ] || mkdir "$comb_dir"
comb_d="$(readlink -fnq "$comb_dir")"
if [ -z $comb_d ]
then
    echo "ERROR: Could not find or create intermediate directory $comb_dir" >&2;
    exit 1;
fi

echo "
#########################
## ORIGINAL DATA: $in_d
## FIGURE INPUT DATA: $out_d
#########################
"

#Process input data into intermediate directory
# see '${combine_script}' for details
echo "
####################################
## Processing original input
###################################
"
bash "$combine_script" "$in_d" "$comb_d" "$ld_file_tag_default" "$expr_file_tag_default" "$expr_signif_file_tag_default" "$rmip_file_tag_default"


#Select subset of data for plotting
# see '${filter_script} for details
echo "
####################################
## Filtering data for figure
###################################
"
bash "$filter_script" "$comb_d" "$out_d" "$comb_d/$ld_file_default" "$comb_d/$expr_file_default" "$comb_d/$expr_signif_file_default" "$comb_d/$rmip_file_default" "$gff_file_source_default" "$bedtools_bin_default"
