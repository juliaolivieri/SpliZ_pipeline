#!/bin/bash
#
#SBATCH --job-name=var
#SBATCH --output=/scratch/PI/horence/JuliaO/single_cell/SZS_pipeline2/scripts/job_output/var.%j.out
#SBATCH --error=/scratch/PI/horence/JuliaO/single_cell/SZS_pipeline2/scripts/job_output/var.%j.err
#SBATCH --time=6:00:00
##SBATCH --qos=normal
#SBATCH -p owners,horence,quake
##SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=40G
date

#DATANAME="TSP2_10x_rerun_with_postprocessing_3prime_cellann"
#DATANAME="TSP1_10x_with_postprocessing_nopanc_cellann"
#DATANAME="TSP2_SS2_RUN1_RUN2_cellann"
DATANAME="TS_pilot_smartseq_with_postprocessing_nopanc_cellann"
a="python3.6 -u /scratch/PI/horence/JuliaO/single_cell/SZS_pipeline2/scripts/variance_adjusted_permutations.py --group_col ontology --suffix _S_0.1_z_0.0_b_5 --dataname ${DATANAME} --num_perms 100 --sub_col tissue"
echo $a 
eval $a
date
