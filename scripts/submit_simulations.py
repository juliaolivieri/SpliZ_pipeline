import glob
import os
import subprocess
import sys
import time
import argparse

def submit_job(file_name):
  """Submit sbatch job to cluster"""
  status, job_num = subprocess.getstatusoutput("sbatch {}".format(file_name))
  if status == 0:
    print("{} ({})".format(job_num, file_name))
    return job_num.split()[-1]
  else:
    print("Error submitting job {} {} {}".format(status, job_num, file_name))

def sbatch_file(file_name,out_path, job_name, time, mem, command, dep="", dep_type = "afterok"):
  """Write sbatch script given parameters"""
  job_file = open(file_name, "w")
  job_file.write("#!/bin/bash\n#\n")
  job_file.write("#SBATCH --job-name=" + job_name + "\n")
  job_file.write("#SBATCH --output={}log_files/{}.%j.out\n".format(out_path, job_name))
  job_file.write("#SBATCH --error={}log_files/{}.%j.err\n".format(out_path, job_name))
  job_file.write("#SBATCH --time={}\n".format(time))
  job_file.write("#SBATCH -p horence,quake\n")
  job_file.write("#SBATCH --nodes=1\n")
  job_file.write("#SBATCH --mem={}\n".format(mem)) 
  if dep != "":
    job_file.write("#SBATCH --dependency={}:{}\n".format(dep_type,dep))
    job_file.write("#SBATCH --kill-on-invalid-dep=yes\n")
  job_file.write("date\n")
  job_file.write(command + "\n")
  job_file.write("date\n")
  job_file.close()

def sig_genes(dataname,z_col,outpath,num):
  command ="python3 -u /scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/significant_genes.py --dataname {} --z_col {} --pinning_S 0.1 --pinning_z 0.0 --lower_bound 5".format(dataname,z_col)
  sbatch_file("{}rewritten_subscripts/run_sig.sh".format(outpath), outpath, "_".join(["sig",num,dataname]),"10:00","60Gb",command)
  return submit_job("{}rewritten_subscripts/run_sig.sh".format(outpath))

def fdr_genes(dataname,z_col,depend,outpath,num):
  command = "python3 -u /scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/final_FDRs_mz.py --dataname {} --suffix _S_0.1_z_0.0_b_5 --z_col {} --nofdr --all_datanames x y".format(dataname,z_col)
  sbatch_file("{}rewritten_subscripts/run_fdr.sh".format(outpath), outpath, "_".join(["fdr",num,dataname]),"10:00","60Gb",command,dep=depend)
  return submit_job("{}rewritten_subscripts/run_fdr.sh".format(outpath))



def box(dataname, letter, params, organism, subfolder, cell_lim, outpath):
  command = "python3 -u /scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/donor_boxplots.py --dataname {} --letter {} --params {} --organism {} --subfolder {} --cell_lim {}".format(dataname, letter, params, organism, subfolder, cell_lim)
  sbatch_file("{}rewritten_subscripts/run_box.sh".format(outpath), outpath, "_".join(["box",letter,dataname]),"24:00:00","80Gb",command)
  return submit_job("{}rewritten_subscripts/run_box.sh".format(outpath))

def dot(params,datanames,organism,subfolder,cell_lim,FDR_col,outpath):
  command = "python3 -u /scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/dotplots_combdataset.py --params {} --datanames {} --organism {} --subfolder {} --cell_lim {} --FDR_col {}".format(params,datanames,organism,subfolder,cell_lim,FDR_col)
  sbatch_file("{}rewritten_subscripts/run_dot.sh".format(outpath), outpath, "dot" + datanames.replace(" ","_"),"24:00:00","80Gb",command)
  return submit_job("{}rewritten_subscripts/run_dot.sh".format(outpath))

def sim(num_trials, num_perms, max_depth, psi, regime,num_cells,outpath):
  command = "python3 -u /scratch/PI/horence/JuliaO/single_cell/SZS_pipeline3/scripts/power_simulation.py --num_trials {} --num_perms {} --max_depth {} --psi {} --regime {} --num_cells {}".format(num_trials, num_perms, max_depth, psi, regime,num_cells)
  sbatch_file("{}rewritten_subscripts/run_sim.sh".format(outpath), outpath, "sim_{}_{}_{}_{}_{}".format(num_trials, num_perms, max_depth, psi, regime),"60:00:00","40Gb",command)
  return submit_job("{}rewritten_subscripts/run_sim.sh".format(outpath))

def main():

  outpath = "/scratch/PI/horence/JuliaO/single_cell/SZS_pipeline3/scripts/"
#  regimes = ["double","double_switch","one_donor","one_donor_fair","casset","exon_skip_twoends"]
#  regimes = ["CD47_type_set1_const","CD47_type_set2_const","CD47_type_set3_const","exon_skip"]
#  regimes = ["casset1","casset2"]
#  regimes = ["CD47_type_set1_const"]
#  regimes = ["double","double_switch"]
#  regimes = ["exon_skip_twoends"]
#  regimes = ["CD47_type_set1"]
  regimes = ["CD47_type_set1","casset1","exon_skip"]
  num_trials = [1000]
  num_perms = [100]
  max_depths = [20]
  num_cells = [20]
  psis = {"double2" : [0.55, 0.6],"double" : [0.55, 0.6,0.65,0.7],"double_switch" : [0.55, 0.6,0.65,0.7],"one_donor" : [0.2],"one_donor_fair" : [0.2], "casset1" : [.65],"casset2" : [.65], "exon_skip_twoends" : [0.55, 0.6,0.65,0.7], "CD47_type" : [.1,.2,.3],"CD47_type_set1" : [0.5],"CD47_type_set2" : [0.5],"CD47_type_set3" : [0.5],"CD47_type_set1_const" : [0],"CD47_type_set2_const" : [0],"CD47_type_set3_const" : [0],"exon_skip":[0.55,0.6,0.65]}
  for regime in regimes:
    for num_trial in num_trials:
      for num_perm in num_perms:
        for max_depth in max_depths:
          for psi in psis[regime]:
            for num_cell in num_cells:
              sim(num_trial, num_perm, max_depth, psi, regime,num_cell, outpath)
              sim(num_trial, num_perm, max_depth, psi, regime + "_null", num_cell,outpath)

main()
