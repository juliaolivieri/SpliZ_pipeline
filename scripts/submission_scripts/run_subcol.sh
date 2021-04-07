#!/bin/bash
#
#SBATCH --job-name=subcol
#SBATCH --output=/scratch/PI/horence/JuliaO/single_cell/SZS_pipeline2/scripts/job_output/subcol.%j.out
#SBATCH --error=/scratch/PI/horence/JuliaO/single_cell/SZS_pipeline2/scripts/job_output/subcol.%j.err
#SBATCH --time=1:00:00
##SBATCH --qos=normal
#SBATCH -p owners,horence,quake
##SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=40G
date

#DATANAME="TSP2_10x_rerun_with_postprocessing_3prime_cellann"
#DATANAME="TSP1_10x_with_postprocessing_nopanc_cellann"
#DATANAME="TSP2_SS2_RUN1_RUN2_cellann"
#DATANAME="HLCA4_P2_10x_with_postprocessing_lung"
#DATANAME="HLCA_smartseq_P3_with_postprocessing_shared"
#DATANAME="HLCA4_P3_10x_with_postprocessing_lung_lungimmuneMacrophage_10"
DATANAME="HLCA4_P2_10x_with_postprocessing_lung_shuffle"
a="python3.6 -u /scratch/PI/horence/JuliaO/single_cell/SZS_pipeline2/scripts/SVD_subcol.py --dataname ${DATANAME} --pinning_S 0.1 --pinning_z 0.0 --lower_bound 5"
echo $a 
eval $a
date
