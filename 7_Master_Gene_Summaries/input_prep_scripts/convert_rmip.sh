#!/bin/sh

#use script.sh input_dir output_directory

in_dir=$1;
out_dir=$2;
rmip_file_tag=$3

if [ -z $rmip_file_tag ] 
then
    $rmip_file_tag='_RMIP'
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
for f in ${in_d}/*${rmip_file_tag}*
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



# files already contain QTL id as last column, need to be concatenated
rmip_out_f="${out_d}/allTraits_RMIP.tsv";

echo -e "\n## COMBINED OUTPUT FILE: $rmip_out_f";

if [ -f "$rmip_f" ] 
then
    echo "WARNING: Output file $rmip_f already exists, appending data to it."
else

    head -qn1 ${in_d}/*${rmip_file_tag}* | \
	sort -u | \
	head -qn1 | \
	sed -re 's/Common\.SI\.ID/QTL_No/' > "$rmip_out_f"
fi
    
for rmip_f in ${in_d}/*${rmip_file_tag}*
do

    echo "# INPUT FILE: $rmip_f"
    tail -n+2 "$rmip_f" >> "$rmip_out_f"
done

    
