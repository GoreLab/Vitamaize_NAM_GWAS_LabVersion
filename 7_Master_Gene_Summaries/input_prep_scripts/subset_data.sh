#!/bin/sh

ld_file_default='region_snps-ld_peak_snp.tsv'
expr_file_default='ExprEffEst_Correl-allQTL.tsv'
rmip_file_default='allTraits_RMIP.tsv'

if [ ${#@} -le 0 ]
then
    echo "Usage: bash $0 input_directory output_directory filter_file [ld_file] [expr_file] [rmip_file]
REQUIRED:
  input_directory: directory containing combined expression file, ld file, and rmip file
  output_directory: directory where to save the subset used for plotting
  filter_file: tab-delimited file containing chromosome, position, and associated QTL ID for the peak markers to use in plotting
OPTIONAL:
  ld_file: filename for ld file
               default = ${ld_file_default}
  expr_file: filename for expression correlation file
               default = ${expr_file_default}
  rmip_file: filename for RMIP file
               default = ${rmip_file_default}
";
    exit;
fi

in_dir=$1; 
out_dir=$2;
filter_file=$3;
ld_file=$4;
expr_file=$5;
rmip_file=$6;

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

filter_f="$(readlink -fnq "$filter_file")"
if [ -z $filter_f ]
then
    echo "ERROR: Could not find filter file $filter_file" >&2;
    exit 1;
fi

if [ -z "$ld_file" ]
then
    ld_file="$in_d/$ld_file_default"
fi
ld_f="$(readlink -fnq "$ld_file")"
if [ -z $ld_f ]
then
    echo "ERROR: Could not find input file $ld_file" >&2;
    exit 1;
fi

if [ -z "$expr_file" ]
then
    expr_file="$in_d/$expr_file_default"
fi
expr_f="$(readlink -fnq "$expr_file")"
if [ -z $expr_f ]
then
    echo "ERROR: Could not find input file $expr_file" >&2;
    exit 1;
fi


if [ -z "$rmip_file" ]
then
    rmip_file="$in_d/$rmip_file_default"
fi
rmip_f="$(readlink -fnq "$rmip_file")"
if [ -z $rmip_f ]
then
    echo "ERROR: Could not find input file $rmip_file" >&2;
    exit 1;
fi


#Expected format of filter file is: "QTL_ID\tChr\tPos"

#Process ld file
echo "Filtering $ld_f"
filtered_ld_file=$(basename "$ld_f")
filtered_ld_file="${filtered_ld_file%.*}-filtered.tsv"
(head -n1 "$ld_f" && \
    cut -f2,3 -d$'\t' "$filter_f" | \
    sed -re 's/\t/./' | \
    sort -ub | \
    join -j1 - <(tail -n+2 "$ld_f" | awk 'BEGIN{FS="\t"}{print $5"."$6"\t"$0}' | sort -k1b,1) -t$'\t' | \
    cut -f2- ) > "$out_d/$filtered_ld_file"


#Process expression file
echo "Filtering $expr_f"
filtered_expr_file=$(basename "$expr_f")
filtered_expr_file="${filtered_expr_file%.*}-filtered.tsv"
(head -n1 "$expr_f" && \
    cut -f1 -d$'\t' "$filter_f" | \
    sort -ub | \
    join -j1 - <(tail -n+2 "$expr_f" | sort -k1b,1 ) -t$'\t' ) > "$out_d/$filtered_expr_file"

#Process RMIP file
echo "Filtering $rmip_f"
filtered_rmip_file=$(basename "$rmip_f")
filtered_rmip_file="${filtered_rmip_file%.*}-filtered.tsv"
(head -n1 "$rmip_f" && \
    cut -f1 -d$'\t' "$filter_f" | \
    sort -ub | \
    join -j1 - <(tail -n+2 "$rmip_f" | sort -k6b,6 ) -t$'\t' -o2.1,2.2,2.3,2.4,2.5,2.6) > "$out_d/$filtered_rmip_file"



