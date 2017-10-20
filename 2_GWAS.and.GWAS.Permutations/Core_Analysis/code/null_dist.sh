#!/bin/bash

#===================================
# SOURCE CONFIG FILE
#++++++++++++++++++++++++++++++++++

#get full (canonical) path to the script location
script_f=$(readlink -enq ${0});
script_d=$(dirname ${script_f});

#source configuration file (required)
#Default: nulldist.cfg in same directory as the script
cfg_f="${script_d}/nulldist.cfg"
if cfg=$(readlink -enq ${cfg_f})
then 
    source ${cfg}
else 
    echo "ERROR:  Configuration file ${cfg_f} was not found!" >&2
    exit 1
fi

#===================================
# FUNCTION DEFINITIONS
#++++++++++++++++++++++++++++++++++

#print instructions for script use when called
function usage () { 
    echo "
Script to generate null distribution for residuals for GWAS analysis using local copy of tassel.
For defaults and other runtime options see configuration file:
${cfg}

REQUIRED PARAMETERS:
  -t <name>
     Name of trait to be analyzed
  -r <filename>
     File containing the residuals for this analysis

OPTIONAL PARAMETERS:
  -i <directory> 
     Directory containing the input file for this analysis
     Default: ${in_d}
  -o <directory>
     Outputs will be saved here in a sub-directory based on trait name
     Default: ${out_d}
  -n <numeric>
     Number of permutations
     Default: ${max_perm}
  -c <numeric>
     Total number of chromosomes to process
     Default: ${nr_chr}
  -s <numeric>
     Starting chromosome
     Default: ${start_chr}
  -m <numeric>
     Map resolution to use in centimorgans)
     Default: ${res}
  -x <directory>
     Directory where the tassel libraries and .jar file are located
     Default: ${tassel_d}
  -p <number>
     Number of concurent permutations to run for each available processor
     Default: ${perm_per_proc}
  -a <number>
     Maximum ram (Gb) available to any child process
     Default: ${max_ram}
  -h Print this help and exit
  -V Verbose mode -- output all parameters before running
  -k use existing output directory and keep existing permutation files
"; 
}

#clean up any child processes on exit

function finish {

    if pgrep -d'|' -P $$ 
    then
        #send SIGTERM to all child processes
        echo -e "Cleaning up child processes:\nSending SIGTERM to child processes of $$\n"
        pkill -15 -P $$
    fi

    #if any still exist wait 10 seconds for cleanup
    if  pgrep -d'|' -P $$   
    then
        sleep 10
    fi

    #then send SIGKILL
    if pgrep -d'|' -P $$ 
    then
        echo -e "Child processes still present, sending SIGKILL"
        pkill -9 -P $$
    fi
}

#set up EXIT trap for cleanup in case of premature exit
trap finish EXIT


#==================================
# PROCESS CONFIG AND OPTIONS
#++++++++++++++++++++++++++++++++++

#print usage if no options were given
if [ -z "$1" ]
then 
    usage
    exit 1
fi

#parse command line options
#NOTE: these will overwrite any options specified in the config file

OPTIND='1'
options=':t:r:i:o:n:c:s:m:x:p:a:hVk'
while getopts $options option
do
    case $option in
        t  ) trait=$OPTARG;;
        r  ) residuals_f=$OPTARG;;
        i  ) in_d=$OPTARG;;
        o  ) out_d=$OPTARG;;
        n  ) max_perm=$OPTARG;;
        c  ) nr_chr=$OPTARG;;
        s  ) start_chr=$OPTARG;;
        m  ) res=$OPTARG;;
	x  ) tassel_d=$OPTARG;;
	k  ) keep='TRUE';;
        p  ) perm_per_proc=$OPTARG;;
        a  ) max_ram=$OPTARG;;
        V  ) verbose='TRUE';;
        h  ) usage; exit;;
        \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done

if [ -z "${trait}" ]
then
    echo "ERROR: Trait name must be specified with -t" >&2
    exit 1
fi

if [ -z "${residuals_f}" ]
then
    echo "ERROR: Residuals file must be specified with -r" >&2
    exit 1
fi

if ! res_f_path=$(readlink -enq ${residuals_f})
then
    echo "ERROR: Cannot find residuals file $residuals_f" >&2
    exit 1
else
    residuals_f=${res_f_path}
fi

# Set the max memory to use for child processes using soft caps in ulimit
#-------------------------------------------------------------------------

ram_request=$((${max_ram}*1048576)) #amount requested by user
local_max_ram=$[ $(free|awk '/^Mem:/{print $2}') - 512000 ] #amount of physical RAM on the system minus 500MB

# If the RAM requested by user exceeds (physical RAM - 500MB), lower it to that value
# We want to avoid splilling over into swap, as that will slow things down tremendously and 
# likely lock up the system for a while.

if [ ${ram_request} -gt ${local_max_ram} ]
 then ram_max=${local_max_ram}
 else ram_max=${ram_request}
fi

# Check to see if our allocation exceeds the system hard limit on resources
v_ulimit_H=$(ulimit -H -v) #hard limit set by the system
if [ ${v_ulimit_H} != 'unlimited' ] && [ ${ram_max} -gt ${v_ulimit_H} ] 
  then ram_max=${v_ulimit_H}
fi

ulimit -S -v ${ram_max}


#Check java heap allocation vs nr of concurrent jobs to spawn
# reduce maxheap to totalmem/nrconcurrent

batchnr=$[${perm_per_proc} * ${nr_proc}]
ram_max_gb=$(printf "%0.0f" $( echo "${ram_max} / 1048576" | bc -l ))
maxheap=$(printf "%0.0f" $( echo "${ram_max_gb} / ${batchnr}" | bc -l))

if [ "$(echo ${maxheap} '<=' 0 | bc -l)" -eq 1 ]
then
    echo "ERROR: Not enough memory! Too many concurrent jobs (${batchnr}) for the ${ram_max_gb}GB max RAM available" >&2
    exit 1
fi

if [ "$(echo ${maxheap} '<' ${Xmx} | bc -l)" -eq 1 ]
then 
    echo "
-----------
 WARNING! 
-----------
         Cannot assign max heap of ${Xmx}GB to Java when running ${batchnr} concurrent jobs with ${ram_max_gb}GB max RAM
         Assigning ${maxheap}GB to max java heap
"
    Xmx=${maxheap}
fi


if [ "$(echo ${maxheap} '<' ${Xms} | bc -l)" -eq 1 ]
then
    echo "
-----------
 WARNING! 
-----------
         Cannot assign a minimum heap of ${Xms}GB to Java when running ${batchnr} concurrent jobs with ${ram_max_gb}GB max RAM
         Assigning ${maxheap}GB to min and max java heap
"
    Xms=${maxheap}
    Xmx=${maxheap}
fi


#set file templates for the chosen resolution
map_chr_t=$(echo -n ${map_t} | sed -re "s/\#\#res\#\#/${res}/g")
founders_chr_t=$(echo -n ${founders_t} | sed -re "s/\#\#res\#\#/${res}/g")
rils_chr_t=$(echo -n ${rils_t} | sed -re "s/\#\#res\#\#/${res}/g")


[ "${verbose}" = 'TRUE' ] && echo "
Starting permutation runs with the following parameters:
  -t ${trait}
  -r ${residuals_f}
  -i ${in_d}
  -o ${out_d}
  -n ${max_perm}
  -c ${nr_chr}
  -s ${start_chr}
  -m ${res}
  -x ${tassel_d}
  -p ${perm_per_proc}
  -a ${max_ram}
  -k ${keep}

System parameters:
  Available processors: ${nr_proc}
  Max RAM available: ${ram_max}
  Java starting heap: ${Xms}g
  Java max heap: ${Xmx}g
  Tassel java libs: ${tassel_d}/${lib_d}
  Tassel jar file: ${tassel_d}/${tassel_jar}

Tassel parameters:
  enterlimit: ${enterlimit}
  iterations: ${iterations}
  maxsnps: ${maxsnps}

Data parameters:
  Map resolution: ${res}

File templates (see config file):
  Map file template: ${map_t}
  Founders file template: ${founders_t}
  RILs file template: ${rils_t}
  Permutation files template: ${perm_t}
" 

#==================================
# MAIN SCRIPT LOGIC
#++++++++++++++++++++++++++++++++++
#
# We're done setting up runtime envirionment and variables
# Next we'll start jobs running tassel
# We assume tassel is not multi-core optimized, so we'll start 
# one job for each processor


#verify tassel location and set java libraries
if ! tassel_d_abs=$(readlink -enq ${tassel_d})
then 
    echo "ERROR: Could not find tassel directory ${tassel_d}" >&2
    exit 1
fi

if ! tassel_lib_d=$(readlink -enq "${tassel_d}/${lib_d}")
then
    echo "ERROR: Could not find tassel library directory ${tassel_d}/${lib_d}" >&2
    exit 1
fi

if ! tassel_exe=$(readlink -enq "${tassel_d}/${tassel_jar}")
then
    echo "ERROR: Could not find tassel jar file ${tassel_d}/${tassel_jar}" >&2
    exit 1
fi

classpath="${tassel_lib_d}/*:${tassel_exe}"

#If output directory exists, make a new one by appending unix time to trait name
#unless -k (keep files) option was passed
if ! trait_out_d=$(readlink -fnq "${out_d}/${trait}_nulldist")
then
    echo "ERROR: Could not find path for output trait directories: ${out_d}" >&2
    exit 1
fi

#If output directory exists, make a new one by appending unix time to trait name
if [ -e ${trait_out_d} ]
then
  if [ ! "$keep" = 'TRUE' ]
  then 
    trait_out_d=${trait_out_d}-$(date +%s)
    if  ! mkdir ${trait_out_d}
      then
      echo "ERROR: Could not make trait output directory: ${trait_out_d}" >&2
      exit 1	
    fi
  fi
else
  if  ! mkdir ${trait_out_d}
    then
    echo "ERROR: Could not make trait output directory: ${trait_out_d}" >&2
    exit 1	
  fi
fi   

echo "Saving output to ${trait_out_d}"


#Create all the permutation files
ncol=$( head -n1 ${residuals_f} | awk -F"\t" '{print NF}' )
cols=$( head -n1 ${residuals_f} )
rows=$( cut -f1 ${residuals_f} | tail -n+2 )

echo "Creating ${max_perm} permutation files"

for i in $( seq 1 ${max_perm} )
do 
    permout_f=$( echo -n ${perm_t} | sed -re "s/\#\#perm\#\#/${i}/g" )
    if ! permout=$(readlink -fnq ${trait_out_d}/${permout_f})
    then
	echo "ERROR: Could not find location for permutation files output ${trait_out_d}." >&2
	exit 1
    fi

    col=2
    cmd='paste'
    while [ ${col} -le ${ncol} ]
    do 
	cmd="${cmd} <(cut -f${col} ${residuals_f} | tail -n+2 | shuf)"
	(( col++ ))
    done
    (echo -e "${cols}" && paste <( echo -e "${rows}") <(eval "${cmd}") ) > $permout
done


#determine the number of concurrent permutations to start at once
max_conc_p=$[${perm_per_proc} * ${nr_proc}]

echo "Running batches of ${max_conc_p} jobs"

#Do permutations one chromosome at a time
chr=${start_chr}
end_chr=$[${start_chr} + ${nr_chr} - 1] #1-based chromosome numbering
while [ "${chr}" -le "${end_chr}" ]
do

    #Set up chromosome files
    map_chr=$(echo -n ${map_chr_t} | sed -re "s/\#\#chr\#\#/${chr}/g")
    founders_chr=$(echo -n ${founders_chr_t} | sed -re "s/\#\#chr\#\#/${chr}/g")
    rils_chr=$(echo -n ${rils_chr_t} | sed -re "s/\#\#chr\#\#/${chr}/g")

    if ! map=$(readlink -enq ${in_d}/${map_chr})
    then
	echo "ERROR: Could not find map file ${in_d}/${map_chr} for chr ${chr}; Trying next chr." >&2
	continue
    fi
    
    if ! founders=$(readlink -enq ${in_d}/${founders_chr})
    then 
	echo "ERROR: Could not find founders file ${in_d}/${founders_chr} for chr ${chr}; Trying next chr." >&2
	continue
    fi
    
    if ! rils=$(readlink -enq ${in_d}/${rils_chr})
    then
	echo "ERROR: Could not find rils file ${in_d}/${rils_chr} for chr ${chr}; Trying next chr." >&2
	continue
    fi

#Set permutation file name base
    model_out_f_base="${trait_out_d}/${trait}_model_chr${chr}_perm"
    steps_out_f_base="${trait_out_d}/${trait}_steps_chr${chr}_perm"
    log_out_f_base="${trait_out_d}/${trait}_log_chr${chr}_perm"

    for p in $( seq 1 ${max_perm} )
    do 
	permin_f=$( echo -n ${perm_t} | sed -re "s/\#\#perm\#\#/${i}/g" )
    if ! permin=$(readlink -fnq ${trait_out_d}/${permin_f})
    then
	echo "ERROR: Could not find location for permutation files input ${trait_out_d}." >&2
	exit 1
    fi
	#set output model file
	if ! model_out_f=$(readlink -fnq "${model_out_f_base}${p}.txt")
	then
	    echo "ERROR: could not set model output file ${model_out_f_base}${p}.txt: cannot resolve path" >&2
	    exit 1
	fi

	#set output steps file
	if ! steps_out_f=$(readlink -fnq "${steps_out_f_base}${p}.txt")
	then
	    echo "ERROR: could not set steps output file ${steps_out_f_base}${p}.txt: cannot resolve path" >&2
	    exit 1
	fi

	#set output log file for each java run
	if ! log_out_f=$(readlink -fnq "${log_out_f_base}${p}.txt")
	then
	    echo "ERROR: could not set log output file ${log_out_f_base}${p}.txt: cannot resolve path" >&2
	    exit 1
	fi


	#=====================================
	#EXECUTION CODE FROM ORIGINAL SCRIPTS
	#-------------------------------------
	# "perl /state/partition1/tassel4-standalone/run_pipeline.pl -Xms500m -Xmx1g -fork1 -NamGwasPlugin -map map.txt -trait trait.txt -rils rils.txt -founders founders.txt -model model.txt -steps steps.txt -enterlimit $enterlimit -iterations $iterations -maxsnps $maxsnps -endPlugin -runfork1"
	#====================================

	java -classpath ${classpath} -Xms${Xms}g -Xmx${Xmx}g net.maizegenetics.pipeline.TasselPipeline -fork1 -NamGwasPlugin -map ${map} -trait ${permin} -rils ${rils} -founders ${founders} -model ${model_out_f} -steps ${steps_out_f} -enterlimit $enterlimit -iterations ${iterations} -maxsnps ${maxsnps} -endPlugin -runfork1 &>${log_out_f} &
	#Line immediately below commented out by CHD Apr 6 2015; uses original residuals file rather than permuted residuals.
	#java -classpath ${classpath} -Xms${Xms}g -Xmx${Xmx}g net.maizegenetics.pipeline.TasselPipeline -fork1 -NamGwasPlugin -map ${map} -trait ${residuals_f} -rils ${rils} -founders ${founders} -model ${model_out_f} -steps ${steps_out_f} -enterlimit $enterlimit -iterations ${iterations} -maxsnps ${maxsnps} -endPlugin -runfork1 &>${log_out_f} &

	#If we reached max concurrent jobs, wait for them to finish
	if [ $(( ${p} % ${max_conc_p} )) -eq 0 ]
	then 
	    wait
	fi

	(( p++ ))

    done

# wait for all spawned jobs to finish before going to the next chromosome
    wait

#concatenate model results into one file/chromosome
    if ! chr_model_out_f=$(readlink -fnq "${trait_out_d}/${trait}_model_chr${chr}_${max_perm}perm_results.txt")
    then
	echo "ERROR: could not set chromosome model output file ${chr_model_out_f}: cannot resolve path" >&2
	exit 1
    fi
    echo -n '' > ${chr_model_out_f}
    cat ${model_out_f_base}* >> ${chr_model_out_f}  

    ((chr++))
done

#run cleanup after everything is done
trap - EXIT
finish
exit 0
