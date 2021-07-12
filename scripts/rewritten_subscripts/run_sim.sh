#!/bin/bash
#
#SBATCH --job-name=sim_1000_100_20_0.65_exon_skip_null
#SBATCH --output=/scratch/PI/horence/JuliaO/single_cell/SZS_pipeline3/scripts/log_files/sim_1000_100_20_0.65_exon_skip_null.%j.out
#SBATCH --error=/scratch/PI/horence/JuliaO/single_cell/SZS_pipeline3/scripts/log_files/sim_1000_100_20_0.65_exon_skip_null.%j.err
#SBATCH --time=60:00:00
#SBATCH -p horence,quake
#SBATCH --nodes=1
#SBATCH --mem=40Gb
date
python3 -u /scratch/PI/horence/JuliaO/single_cell/SZS_pipeline3/scripts/power_simulation.py --num_trials 1000 --num_perms 100 --max_depth 20 --psi 0.65 --regime exon_skip_null --num_cells 20
date
