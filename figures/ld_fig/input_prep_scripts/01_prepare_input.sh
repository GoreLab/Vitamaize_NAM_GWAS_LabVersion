#!/bin/sh

ld_file_tag_default='.target_ld.'
expr_file_tag_default='Expr.EffEst_Correl_QTL'
expr_signif_file_tag_default='Signif_of_Expr.EffEst.Correls'
rmip_file_tag_default='_RMIPgt4_with_QTLnumber_forFig2'

if [ ${#@} -le 0 ]
then
    echo "Usage: bash $0 input_directory output_directory [ld_file_tag] [expr_file_tag] [expr_signif_file_tag] [rmip_file_tag]
REQUIRED:
  input_directory: directory containing expression data, ld data, and gwas data subdirectories
  output_directory: directory where to save transformed data that will be used as input to plotting
OPTIONAL:
  ld_file_tag: unique filename pattern for ld files
               default = ${ld_file_tag_default}
  expr_file_tag: unique filename pattern for expression correlation r-value files
               default = ${expr_file_tag_default}
  expr_signif_file_tag: unique filename pattern for expression correlation significance
               default = ${expr_signif_file_tag_default}
  rmip_file_tag: unique filename pattern for RMIP file(s)
               default = ${rmip_file_tag_default}
";
    exit;
fi

in_dir=$1; #directory containing expression data, ld data, and gwas data subdirectories
out_dir=$2; #location of transformed data to use as input to plotting
ld_file_tag=$3;
expr_file_tag=$4;
expr_signif_file_tag=$5;
rmip_file_tag=$6;

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

if [ -z "$ld_file_tag" ]
then
    ld_file_tag="$ld_file_tag_default"
fi

if [ -z "$expr_file_tag" ]
then
    expr_file_tag="$expr_file_tag_default"
fi

if [ -z "$expr_signif_file_tag" ]
then
    expr_signif_file_tag="$expr_signif_file_tag_default"
fi

if [ -z "$rmip_file_tag" ]
then
    rmip_file_tag="$rmip_file_tag_default"
fi


# build a list of all subdirectories with *.target_ld.* files, and process each subdirectory sequentially
for ld_dir in $( for f in ${in_d}/*/*${ld_file_tag}* ; do echo $(dirname "$f") ; done | sort -u );
do
    echo -e "\n\n########################"
    echo "## Processing ${ld_dir}"
#append peak location info as last columns and combine all files
    bash convert_ldsnps.sh "$ld_dir" "$out_d" "$ld_file_tag"
done

# build a list of all subdirectories with *_RMIP* files, and process each subdirectory sequentially
for rmip_dir in $( for f in ${in_d}/*/*${rmip_file_tag}* ; do echo $(dirname "$f") ; done | sort -u );
do
    echo -e "\n\n########################"
    echo "## Processing ${rmip_dir}"
    # file already contains QTL id as last column
    bash convert_rmip.sh "$rmip_dir" "$out_d" "$rmip_file_tag"
done


# build a list of all subdirectories with *Expr.EffEst_Correl* files, and process each subdirectory sequentially
for expr_dir in $( for f in ${in_d}/*/*${expr_file_tag}* ; do echo $(dirname "$f") ; done | sort -u );
do
    echo -e "\n\n########################"
    echo "## Processing ${expr_dir}"
#insert QTL id as first column and combine all files
#this data format will be melted in R later
    bash convert_expr.sh "$expr_dir" "$out_d" "$expr_file_tag"
done


# build a list of all subdirectories with significance values for expr correlations, and process each subdirectory sequentially
for expr_signif_dir in $( for f in ${in_d}/*/*${expr_signif_file_tag}* ; do echo $(dirname "$f") ; done | sort -u );
do
    echo -e "\n\n########################"
    echo "## Processing ${expr_signif_dir}"
#remove rows with no data and move to output
#this data format will be melted in R later
    bash convert_expr_sig.sh "$expr_signif_dir" "$out_d" "$expr_signif_file_tag"
done

