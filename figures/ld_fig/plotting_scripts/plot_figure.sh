#!/bin/sh

#===================================
# DEFAULT VALUES
#+++++++++++++++++++++++++++++++++++

# I/O defaults
output_file_name_default='Fig2'
output_file_type_default='pdf'
ld_file_name_default='fig_data-ld_from_peaks.tsv';
expr_corr_file_name_default='fig_data-gene_expr_correl.tsv';
expr_corr_sig_file_name_default='fig_data-gene_expr_correl_sig.tsv'
rmip_file_name_default='fig_data-gwas_rmip.tsv';
gene_coord_file_name_default='fig_data-candidate_gene_coord.bed';
key_file_name_default='fig_data-combined_filter_list.tsv';

# runtime parameters
search_space_default='100000'
bin_size_default='2000'
pval_cutoff_default='0.05'
ribbon_type_default='flat'
label_type_default='arrow'




#===================================
# FUNCTION DEFINITIONS
#+++++++++++++++++++++++++++++++++++

#print instructions for script use when called
usage () {

    echo "Usage: bash $0 -i input_directory -o output_directory [ -f output_file_name ] [ -t output_file_type ] [ -l ld_file_name ] [ -e expr_corr_file_name ] [ -s expr_corr_sig_file_name ] [ -r rmip_file_name ] [ -g gene_coord_file_name ] [ -k key_file_name ] [ -d search_space ] [ -b bin_size ] [ -p pval_cutoff ] [ -R ribbon_type ] [ -L label_type ]
REQUIRED:
  input_directory: directory containing data filtered for plotting
  output_directory: directory where to save data figures output
OPTIONAL:
  output_file_name: filename (no path or suffix) for the output file
               default: '${output_file_name_default}'
  output_file_type: type of output file, one of: png, jpeg, tiff, pdf, svg, ps
               default: '${output_file_type_default}'
  ld_file_name: filename for ld input file
               default = '${ld_file_name_default}'
  expr_corr_file_name: filename for expression correlation file
               default = '${expr_corr_file_name_default}'
  expr_corr_sig_file_name: filename for expression correlation significance file
               default = '${expr_corr_sig_file_name_default}'
  rmip_file_name: filename for gwas RMIP file
               default = '${rmip_file_name_default}'
  gene_coord_file_name: filename for candidate genes .bed file
               default = '${gene_coord_file_name_default}'
  key_file_name: filename for the keyfile linking QTL, gene, and peak SNP
               default = '${key_file_name_default}'
  search_space: distance (bp) on either side of peak to show for figure
               default = '${search_space_default}'  
  bin_size: size of bins (bp)
               default = '${bin_size_default}'
  pval_cutoff: cutoff pvalue for significant expression correlation
               default = '${pval_cutoff_default}'
  ribbon_type: switch for graphical display of ld bins; one of 'flat', 'hist'
               default = '${ribbon_type_default}'
  label_type:  switch for labeling of GWAS RMIP values; one of 'arrow', 'text', 'off'
               default = '${label_type_default}'
";
}



if [ ${#@} -le 0 ]
then
    usage
    exit 1;
fi


#parse command line options
OPTIND='1'
options=':i:o:f:t:l:e:s:r:g:k:d:b:p:R:L:Vh'
while getopts $options option
do
    case $option in
	i  ) data_dir=$OPTARG;;
	o  ) fig_dir=$OPTARG;;
	f  ) output_file_name=$OPTARG;;
	t  ) output_file_type=$OPTARG;;
	l  ) ld_file_name=$OPTARG;;
	e  ) expr_file_name=$OPTARG;;
        s  ) expr_sig_file_name=$OPTARG;;
	r  ) rmip_file_name=$OPTARG;;
	g  ) genes_file_name=$OPTARG;;
	k  ) key_file_name=$OPTARG;;
	d  ) search_space=$OPTARG;;
	b  ) bin_size=$OPTARG;;
	p  ) pval_cutoff=$OPTARG;;
	R  ) ribbon_type=$OPTARG;;
        L  ) label_type=$OPTARG;;
        V  ) verbose='TRUE';;
        h  ) usage; exit;;
        \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done




data_d="$(readlink -enq "$data_dir")"
if [ -z "$data_d" ]
then
    echo "ERROR: Could not find input data directory $data_dir" >&2;
    exit 1
fi

fig_d="$(readlink -enq "$fig_dir")"
if [ -z "$fig_d" ]
then
    echo "ERROR: Could not find output data directory $fig_dir" >&2;
    exit 1
fi


echo -e "Reading data from ${data_d}\nSaving output to ${fig_d}"


if [ -z "$output_file_type" ]
then
    output_file_type="$output_file_type_default"
fi
if [ -z "$output_file_name" ]
then
    output_file="$fig_d/${output_file_name_default}.${output_file_type}"
else
    output_file="$fig_d/${output_file_name}.${output_file_type}"
fi
output_f="$(readlink -fnq "$output_file")"
if [ -z "$output_f" ]
then
    echo "ERROR: Could not find path to output file ${output_file}" >&2;
    exit 1
fi




if [ -z "$ld_file_name" ]
then
    ld_file="$data_d/$ld_file_name_default"
else
    ld_file="$data_d/$ld_file_name"
fi
ld_f="$(readlink -enq "$ld_file")"
if [ -z "$ld_f" ]
then
    echo "ERROR: Could not find input file ${ld_file}" >&2;
    exit 1
fi


if [ -z "$expr_file_name" ]
then
    expr_file="$data_d/$expr_corr_file_name_default"
else
    expr_file="$data_d/$expr_file_name"
fi
expr_f="$(readlink -enq "$expr_file")"
if [ -z "$expr_f" ]
then
    echo "ERROR: Could not find input file ${expr_file}" >&2;
    exit 1
fi

if [ -z "$expr_sig_file_name" ]
then
    expr_sig_file="$data_d/$expr_corr_sig_file_name_default"
else
    expr_sig_file="$data_d/$expr_sig_file_name"
fi
expr_sig_f="$(readlink -enq "$expr_sig_file")"
if [ -z "$expr_sig_f" ]
then
    echo "ERROR: Could not find input file ${expr_sig_file}" >&2;
    exit 1
fi


if [ -z "$rmip_file_name" ]
then
    rmip_file="$data_d/$rmip_file_name_default"
else
    rmip_file="$data_d/$rmip_file_name"
fi
rmip_f="$(readlink -enq "$rmip_file")"
if [ -z "$rmip_f" ]
then
    echo "ERROR: Could not find input file ${rmip_file}" >&2;
    exit 1
fi


if [ -z "$genes_file_name" ]
then
    genes_file="$data_d/$gene_coord_file_name_default"
else
    genes_file="$data_d/$genes_file_name"
fi
genes_f="$(readlink -enq "$genes_file")"
if [ -z "$genes_f" ]
then
    echo "ERROR: Could not find input file ${genes_file}" >&2;
    exit 1
fi

if [ -z "$key_file_name" ]
then
    key_file="$data_d/$key_file_name_default"
else
    key_file="$data_d/$key_file_name"
fi
key_f="$(readlink -enq "$key_file")"
if [ -z "$key_f" ]
then
    echo "ERROR: Could not find input file ${key_file}" >&2;
    exit 1
fi

if [ -z "$search_space" ]
then
    search_space=$search_space_default
fi

if [ -z "$bin_size" ]
then
    bin_size=$bin_size_default
fi

if [ -z "$pval_cutoff" ]
then
    pval_cutoff=$pval_cutoff_default
fi

if [ -z "$ribbon_type" ]
then
    ribbon_type=$ribbon_type_default
fi

if [ -z "$label_type" ]
then
    label_type=$label_type_default
fi



export fig_d
export output_f
export output_file_type
export ld_f
export expr_f
export expr_sig_f
export rmip_f
export genes_f
export key_f
export search_space
export bin_size
export pval_cutoff

echo "Plotting and saving figure to ${output_f}"
R --slave --vanilla -f "$(dirname "$(readlink -fnq "$0")")/fig2.R"
