#!/bin/bash

#===================================
# SOURCE CONFIG FILE
#+++++++++++++++++++++++++++++++++++

#get full (canonical) path to the script location
script_f=$(readlink -enq ${0});
script_d=$(dirname ${script_f});

#source configuration file (required)
#Default: iterations.cfg in same directory as the script
cfg_f="${script_d}/iterations.cfg"
if cfg=$(readlink -enq ${cfg_f})
then 
    source ${cfg}
else 
    echo "ERROR:  Configuration file ${cfg_f} was not found!" >&2
    exit 1
fi

#===================================
# FUNCTION DEFINITIONS
#+++++++++++++++++++++++++++++++++++

#print instructions for script use when called
usage () { 
    echo "
Script to run several permutations of GWAS analysis using local copy of tassel.
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
     Number of itterations to run
     Default: ${total_iterations}
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
     Number of concurrent tassel runs
     Default: ${nr_proc}
  -a <number>
     Maximum ram (Gb) available to any child process
     Default: ${max_ram}
  -h Print this help and exit

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
        pkill -15 -P $$
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
options=':t:r:i:o:n:c:s:m:x:p:a:h'
while getopts $options option
do
    case $option in
        t  ) trait=$OPTARG;;
        r  ) residuals_f=$OPTARG;;
        i  ) in_d=$OPTARG;;
        o  ) out_d=$OPTARG;;
        n  ) total_iterations=$OPTARG;;
        c  ) nr_chr=$OPTARG;;
        s  ) start_chr=$OPTARG;;
        m  ) res=$OPTARG;;
	x  ) tassel_d=$OPTARG;;
	p  ) nr_proc=$OPTARG;;
        a  ) max_ram=$OPTARG;;
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

#set file templates for the chosen resolution
map_chr_t=$(echo -n ${map_t} | sed -re "s/\#\#res\#\#/${res}/g")
founders_chr_t=$(echo -n ${founders_t} | sed -re "s/\#\#res\#\#/${res}/g")
rils_chr_t=$(echo -n ${rils_t} | sed -re "s/\#\#res\#\#/${res}/g")


###########################################################
# We're done setting up runtime envirionment and variables
# Next we'll start jobs running tassel
# We assume tassel is not multi-core optimized, so we'll start 
# one job for each processor

echo "
Starting permutation runs with the following parameters:
  -t ${trait}
  -r ${residuals_f}
  -i ${in_d}
  -o ${out_d}
  -n ${total_iterations}
  -c ${nr_chr}
  -s ${start_chr}
  -m ${res}
  -x ${tassel_d}
  -p ${nr_proc}
  -a ${max_ram}

Other variables:
  Java starting heap: ${Xms}g
  Java max heap: ${Xmx}g
  Tassel java libs: ${tassel_d}/${lib_d}
  Tassel jar file: ${tassel_d}/${tassel_jar}
  enterlimit: ${enterlimit}
  Map resolution: ${res}

File templates (see config file):
  Map file template: ${map_t}
  Founders file template: ${founders_t}
  RILs file template: ${rils_t}
" 


#==================================
# MAIN SCRIPT LOGIC
#++++++++++++++++++++++++++++++++++


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


#make output directory for the trait
if ! trait_out_d=$(readlink -fnq "${out_d}/${trait}_iterations")
then
    echo "ERROR: Could not find path for output trait directories: ${out_d}" >&2
    exit 1
fi

#If output directory exists, make a new one by appending unix time to trait name
if [ -e ${trait_out_d} ]
then
    trait_out_d=${trait_out_d}-$(date +%s)
fi

if ! mkdir ${trait_out_d}
then
    echo "ERROR: Could not make trait output directory: ${trait_out_d}" >&2
    exit 1
fi

echo "Saving output to ${trait_out_d}"




#Divide iterations among the number of concurent runs
#If division has remainder, the last run assigned will have less than ${step} itterations
iter_step=$(perl -e 'use POSIX; print ceil($ARGV[0]/$ARGV[1]), qq{\n}' ${total_iterations} ${nr_proc})


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
    
    iter_start='0'
    while [ "${iter_start}" -lt "${total_iterations}" ]
    do
	step=${iter_step}
	iter_end=$[${iter_start} + ${iter_step}]
	if [ "${iter_end}" -gt "${total_iterations}" ]
	then
	    iter_end=${total_iterations}
	    step=$[${iter_end} - ${iter_start}]
	fi

	#set output model file
	if ! model_out_f=$(readlink -fnq "${trait_out_d}/${trait}_model_chr${chr}_iter${iter_start}-${iter_end}.txt")
	then
	    echo "ERROR: could not set model output file ${model_out_f}: cannot resolve path" >&2
	    exit 1
	fi

	#set output steps file
	if ! steps_out_f=$(readlink -fnq "${trait_out_d}/${trait}_steps_chr${chr}_iter${iter_start}-${iter_end}.txt")
	then
	    echo "ERROR: could not set steps output file ${steps_out_f}: cannot resolve path" >&2
	    exit 1
	fi

	#set output log file for each java run
	if ! log_out_f=$(readlink -fnq "${trait_out_d}/${trait}_chr${chr}_iter${iter_start}-${iter_end}-log.txt")
	then
	    echo "ERROR: could not set log output file ${log_out_f}: cannot resolve path" >&2
	    exit 1
	fi

	#=====================================
	#EXECUTION CODE FROM ORIGINAL SCRIPTS
	#-------------------------------------
	#echo "perl ${indir}/tassel4-updated/run_pipeline.pl -Xms500m -Xmx800m -fork1 -NamGwasPlugin -map ${map} -trait ${residuals} -rils ${rils} -founders ${founders} -model ${outdir}/${trait}/${trait}_model_chr${chr}_iter${iter}.txt -steps ${outdir}/${trait}/${trait}_steps_chr${chr}_iter${iter}.txt -enterlimit $enterlimit -iterations $step -start $iter -fullmodel -endPlugin -runfork1 &" >> $script
	#print "java -classpath '$CP' $java_mem_min $java_mem_max net.maizegenetics.pipeline.TasselPipeline @args\n";
	#====================================


	java -classpath ${classpath} -Xms${Xms}g -Xmx${Xmx}g net.maizegenetics.pipeline.TasselPipeline -fork1 -NamGwasPlugin -map ${map} -trait ${residuals_f} -rils ${rils} -founders ${founders} -model ${model_out_f} -steps ${steps_out_f} -enterlimit $enterlimit -iterations ${step} -start ${iter_start} -fullmodel -endPlugin -runfork1 &>${log_out_f} &
	(( iter_start += ${step} ))
	(( iter_start++ ))
    done

# wait for all spawned jobs to finish before going to the next chromosome
    wait

    ((chr++))
done

#run cleanup after everything is done
trap - EXIT
finish
exit 0

