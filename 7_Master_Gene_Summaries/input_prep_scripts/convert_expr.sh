#!/bin/sh

#use script.sh in_directory output_directory

in_dir=$1;
out_dir=$2;
expr_file_tag=$3;

if [ -z $expr_file_tag ]
then
    expr_file_tag='Expr.EffEst'
fi


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


# make sure files are in proper unix text format
for f in ${in_d}/*${expr_file_tag}* 
do 
    f_form=$(file -b "$f") 
    ftype=${f_form%%,*}
    flines=$(echo ${f_form##*,} | grep -woP '(CR|LF|CRLF)')
    if [ ! "${ftype%% *}" == 'ASCII' ]
    then 
	echo "ERROR: $f is not an ASCII text file" >&2
    else 
	if [ -n "$flines" ] && [ "$flines" == 'CRLF' ]
	then 
	    dos2unix $f
	else 
	    if [ -n "$flines" ] && [ "$flines" == 'CR' ]
	    then 
		mac2unix $f
	    fi
	fi
    fi
done

# insert qtl id as the first column and combine all qtl files
expr_out_f="${out_d}/ExprEffEst_Correl-allQTL.tsv";

echo -e "\n## COMBINED OUTPUT FILE: $expr_out_f"

if [ -f "$expr_out_f" ] 
then
    echo "WARNING: Output file $expr_out_f already exists, appending data to it."
else
    
    head -qn1 ${in_d}/*${expr_file_tag}* | \
	sort -u | \
	head -qn1 | \
	sed -re 's/^/QTL_No\tQTL\t/' -e 's/\t([0-9]+)/\t\1DAP/g' > "$expr_out_f"
fi
for expr_f in ${in_d}/*${expr_file_tag}*
do 
    qtlID=$(basename "$expr_f")
    qtlID=${qtlID##*_QTL}
    qtlID=${qtlID%.*}

    echo "# INPUT FILE: $expr_f"
    echo "# QTL: $qtlID"
echo "# QTL NR: ${qtlID%_*}"
    tail -n+2 "$expr_f" | \
	grep -vP '^\s*$' | \
	sed -re "s/^/${qtlID%_*}\t${qtlID}\t/" | \
	sort -k1n,1 -k2V,2 -k3,3 >> "$expr_out_f"
done 


