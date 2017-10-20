#!/bin/bash
#
ulimit -Su $(ulimit -Hu)
#bash null_dist.sh -t PC8 -r ../../data_input/Resid.PC8.SI.01_2015.txt -p 2 -i ../../data_input/ -o /home/mag87_0001/chd45/Documents/permutation_scripts-v3_1-20150402/analysis_output/ -a 500 #run completed on CBSU server
#sleep 2m
bash null_dist.sh -t totalTocochrs -r ../../data_input/Resid.totalTocochrs.SI.01_2015.txt -p 1 -i ../../data_input/ -o ../analysis_output/Complete_Results_BashPart3/ -a 13 -s 3 -c 1 -V -k #chrs 1 & 2 finished on CBSU server, copied partial directory to lab server for resuming; note 4/18 looked for chr 11 because didn't specify # of chrs; results should still be ok, ran to completion for chr 10. 4/22 re-ran chr 3 so that results obtained without appending to model output files
sleep 2m
#bash null_dist.sh -t ACAR_RUV -r ../../data_input/Resid.ACAR_RUV.SI.01_2015.txt -p 1 -i ../../data_input/ -o ../analysis_output/Complete_Results_BashPart3/ -V -a 13 #finished on lab server 4/21 3am
#sleep 2m
#bash null_dist.sh -t BCAR_RUV -r ../../data_input/Resid.BCAR_RUV.SI.01_2015.txt -p 1 -i ../../data_input/ -o ../analysis_output/Complete_Results_BashPart3/ -V -a 13 -s 1 -c 1 -k #CD added s 2 & c9 on 4/21 to finish runs stopped upon power outage, but made new directory (& thus made new perm files) to avoid overwrite. Dir changed to _nulldist/, then re-ran -c 1 -s 1 with KEEP to run for chr1
#sleep 2m