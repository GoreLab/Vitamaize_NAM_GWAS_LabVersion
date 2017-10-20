#!/bin/bash
#
bash run_iterations.sh -t aT -r /share/maize_gwas/data_input/Resid.aT.SI.01.txt -e 2.48e-06 -i /share/maize_gwas/data_input -o /home/chd45/NAM_GWAS/GWAS_Results_Thresholds.by.Trait/aT_iterations/ -p 3 -a 13 -V &
bash run_iterations.sh -t aT3 -r /share/maize_gwas/data_input/Resid.aT3.SI.01.txt -e 2.31e-06 -i /share/maize_gwas/data_input -o /home/chd45/NAM_GWAS/GWAS_Results_Thresholds.by.Trait/aT3_iterations/ -p 3 -a 13 -V &
bash run_iterations.sh -t dT -r /share/maize_gwas/data_input/Resid.dT.SI.01.txt -e 3.50e-06 -i /share/maize_gwas/data_input -o /home/chd45/NAM_GWAS/GWAS_Results_Thresholds.by.Trait/dT_iterations/ -p 3 -a 13 -V &
bash run_iterations.sh -t dT3 -r /share/maize_gwas/data_input/Resid.dT3.SI.01.txt -e 9.22e-07 -i /share/maize_gwas/data_input -o /home/chd45/NAM_GWAS/GWAS_Results_Thresholds.by.Trait/dT3_iterations/ -p 3 -a 13 -V &
bash run_iterations.sh -t gT -r /share/maize_gwas/data_input/Resid.gT.SI.01.txt -e 2.38e-06 -i /share/maize_gwas/data_input -o /home/chd45/NAM_GWAS/GWAS_Results_Thresholds.by.Trait/gT_iterations/ -p 3 -a 13 -V
sleep 2m
bash run_iterations.sh -t gT3 -r /share/maize_gwas/data_input/Resid.gT3.SI.01.txt -e 4.16e-06 -i /share/maize_gwas/data_input -o /home/chd45/NAM_GWAS/GWAS_Results_Thresholds.by.Trait/gT3_iterations/ -p 3 -a 13 -V &
bash run_iterations.sh -t PC8 -r /share/maize_gwas/data_input/Resid.PC8.SI.01.txt -e 3.65e-06 -i /share/maize_gwas/data_input -o /home/chd45/NAM_GWAS/GWAS_Results_Thresholds.by.Trait/PC8_iterations/ -p 3 -a 13 -V &
bash run_iterations.sh -t totalT -r /share/maize_gwas/data_input/Resid.totalT.SI.01.txt -e 6.11e-06 -i /share/maize_gwas/data_input -o /home/chd45/NAM_GWAS/GWAS_Results_Thresholds.by.Trait/totalT_iterations/ -p 3 -a 13 -V &
bash run_iterations.sh -t totalT3 -r /share/maize_gwas/data_input/Resid.totalT3.SI.01.txt -e 3.37e-06 -i /share/maize_gwas/data_input -o /home/chd45/NAM_GWAS/GWAS_Results_Thresholds.by.Trait/totalT3_iterations/ -p 3 -a 13 -V &
bash run_iterations.sh -t totalTocochrs -r /share/maize_gwas/data_input/Resid.totalTocochrs.SI.01.txt -e 3.77e-06 -i /share/maize_gwas/data_input -o /home/chd45/NAM_GWAS/GWAS_Results_Thresholds.by.Trait/totalTocochrs_iterations/ -p 3 -a 13 -V
sleep 2m