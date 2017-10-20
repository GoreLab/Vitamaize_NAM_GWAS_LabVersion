#!/bin/bash
#
#bash null_dist.sh -t aT -r ../../data_input/Resid.aT.SI.01_2015.txt -p 1 -i ../../data_input/ -a 500 #run on CBSU server Sat to Sun 4/11-12/2015 CHD
#sleep 2m
#cp -r -n analysis_output /home/mag87_0001/chd45/
bash null_dist.sh -t dT3 -r ../../data_input/Resid.dT3.SI.01_redo_ExtrValRemoved_newTransAndPerm.txt -p 1 -i ../../data_input/ -o ../analysis_output/Complete_Results_dT3_ExtrValRemoved_newTransAndPerm/ -a 13 -s 1 -c 10 -V -k #run on CBSU server Sun to Mon 4/12-13/2015 CHD, restarted on chr 10 using same perm files (edit made by Dan 4/13). 4/22 re-ran chr 10 on lab server so that results (model output) would not append, rather only 1 line per file.
sleep 2m
#cp -r -n analysis_output /home/mag87_0001/chd45/
#ulimit -Su $(ulimit -Hu)
#bash null_dist.sh -t dT -r ../../data_input/Resid.dT.SI.01_2015.txt -p 2 -i ../../data_input/ -a 500
#sleep 2m
#cp -r -n ../analysis_output /home/mag87_0001/chd45/Documents/permutation_scripts-v3_1-20150402/analysis_output

#bash null_dist.sh -t dT3 -r ../../data_input/Resid.dT3.SI.01_2015.txt -p 2 -i ../../data_input/ -a 500
#sleep 2m
#cp -r -n ../analysis_output /home/mag87_0001/chd45/Documents/permutation_scripts-v3_1-20150402/analysis_output
