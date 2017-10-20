#!/bin/bash
#
bash null_dist.sh -t TOTCAR_RUV -r ../../data_input/Resid.TOTCAR_RUV.SI.01_2015.txt -p 3 -i ../../data_input/
sleep 2m
bash null_dist.sh -t ZEA_RUV -r ../../data_input/Resid.ZEA_RUV.SI.01_2015.txt -p 3 -i ../../data_input/
sleep 2m
bash null_dist.sh -t ZEI_RUV -r ../../data_input/Resid.ZEI_RUV.SI.01_2015.txt -p 3 -i ../../data_input/
sleep 2m

#cp -r -n analysis_output /home/mag87_0001/chd45/ #commented out CHD 4/10 because running part 5 on lab server