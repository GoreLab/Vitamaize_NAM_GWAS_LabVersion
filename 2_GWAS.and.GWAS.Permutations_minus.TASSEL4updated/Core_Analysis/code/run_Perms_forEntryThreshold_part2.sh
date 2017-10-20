#!/bin/bash
#
bash null_dist.sh -t gT -r ../../data_input/Resid.gT.SI.01_2015.txt -p 3 -i ../../data_input/
sleep 2m
bash null_dist.sh -t gT3 -r ../../data_input/Resid.gT3.SI.01_2015.txt -p 3 -i ../../data_input/
sleep 2m
bash null_dist.sh -t totalT ../../data_input/Resid.totalT.SI.01_2015.txt -p 3 -i ../../data_input/
sleep 2m
bash null_dist.sh -t totalT3 -r ../../data_input/Resid.totalT3.SI.01_2015.txt -p 3 -i ../../data_input/
sleep 2m

#cp -r -n analysis_output /home/mag87_0001/chd45/ #this file run on lab server beg. 4/13/2015