#!/bin/sh

#use script.sh input_dir output_directory

in_dir=$1;
out_dir=$2;
ld_file_tag=$3;

if [ -z $ld_file_tag ]
then
    ld_file_tag='.target_ld.'
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
for f in ${in_d}/*${ld_file_tag}*
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


# append peakSNP location as the last column and combine all SNP ld data into one file
linked_snps_f="${out_d}/region_snps-ld_peak_snp.tsv";

echo -e "\n## COMBINED OUTPUT FILE: $linked_snps_f"

if [ -f "$linked_snps_f" ] 
then
    echo "WARNING: Output file $linked_snps_f already exists, appending data to it."
else

    echo -e "Chr\tPos\tAlleles\tRsq\tPeakChr\tPeakPos\tPeakAllele" > "$linked_snps_f";
fi
    
for rsq_f in ${in_d}/*${ld_file_tag}*
do

    peak_info=$(basename "$rsq_f");
    peak_info=${peak_info%%.*};
    peak_info=${peak_info#*+};
    peak_info=$(echo -n "$peak_info" | sed -re 's/(D|[0-9])\+/\1\t/g');

    echo "# INPUT FILE: $rsq_f";
    echo "# PEAK: $peak_info";
    
    tail -n+2 "$rsq_f" | \
	cut -f 1-3,5 | \
	sed -re "s/$/\t${peak_info}/" >> "$linked_snps_f"
done

    
