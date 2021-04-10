datasets = ["HLCA4_P2_10x_with_postprocessing_lung","HLCA4_P3_10x_with_postprocessing_lung"]

num_perms = 100

# filter by SICILIAN v2
ver = "--v2"

#all_groups = {
#              "HLCA4_P2_10x_with_postprocessing_lung" : ["HLCA4_P2_10x_with_postprocessing_lung","HLCA4_P3_10x_with_postprocessing_lung"],
#              "HLCA4_P3_10x_with_postprocessing_lung" : ["HLCA4_P2_10x_with_postprocessing_lung","HLCA4_P3_10x_with_postprocessing_lung"]
#}

light = False
z_col = "scZ"

if light:
  suff = "_light"
else:
  suff = ""

unfilt = False
if unfilt:
  suff += "_unfilt"

corr_methods = ["pearson","spearman"]

verbose = "--verbose"
pins_S = [0.1]
pins_z = [0.0]
bounds = [5]

def get_anova(datasets):
  out = []
  for dataset in datasets:
    for b in bounds:
      for pin_S in pins_S:
        for pin_z in pins_z:
          out.append("scripts/output/anova_zscore/{}_{}_anova_out_S_{}_z_{}_b_{}{}.tsv".format(dataset,z_col,pin_S,pin_z,b,suff))
  return out

def get_FDR(datasets):
  out = []
  for dataset in datasets:
    for b in bounds:
      for pin_S in pins_S:
        for pin_z in pins_z:
          out.append("scripts/output/final_FDRs_mz/{}_FDR_S_{}_z_{}_b_{}{}.tsv".format(dataset,pin_S,pin_z,b,suff))
          out.append("scripts/output/final_FDRs_anova/{}_FDR_S_{}_z_{}_b_{}{}.tsv".format(dataset,pin_S,pin_z,b,suff))
#          pass

  return out

def get_SVD(datasets):
  out = []
  for b in bounds:
    for pin_S in pins_S:
      for pin_z in pins_z:
        for dataset in datasets:
          out.append("scripts/output/rijk_zscore/{}_sym_SVD_normdonor_S_{}_z_{}_b_{}{}.tsv".format(dataset,pin_S,pin_z,b,suff))

  return out

def get_rijk_zscores(datasets):
  out = []
  for b in bounds:
    for pin_S in pins_S:
      for pin_z in pins_z:
        for dataset in datasets:
          out.append("scripts/output/rijk_zscore/{}_sym_S_{}_z_{}_b_{}{}.tsv".format(dataset,pin_S,pin_z,b,suff))

  return out

def get_sig(datasets, bounds):
  out = []
  for b in bounds:
    for pin_S in pins_S:
      for pin_z in pins_z:
        for dataset in datasets:
          out.append("scripts/output/significant_genes/{}-{}_allp_S_{}_z_{}_b_{}{}.tsv".format(dataset,z_col,pin_S,pin_z,b, suff))

  return out

def get_group_anova(wildcards):
  ins = ["scripts/output/anova_zscore/{}_{}_anova_out_S_{}_z_{}_b_{}{}.tsv".format(wildcards.dataset,z_col,wildcards.pinS,wildcards.pinz,wildcards.bound,suff),"scripts/output/anova_zscore/{}_{}_anova_coeff_S_{}_z_{}_b_{}{}.tsv".format(wildcards.dataset,z_col,wildcards.pinS,wildcards.pinz,wildcards.bound,suff)]
  for name in all_groups[wildcards.dataset]:
    ins.append("scripts/output/anova_zscore/{}_{}_anova_out_S_{}_z_{}_b_{}{}.tsv".format(name,z_col,wildcards.pinS,wildcards.pinz,wildcards.bound,suff))
    ins.append("scripts/output/anova_zscore/{}_{}_anova_coeff_S_{}_z_{}_b_{}{}.tsv".format(name,z_col,wildcards.pinS,wildcards.pinz,wildcards.bound,suff))

  return list(set(ins))

def get_group(wildcards):
  ins = ["scripts/output/significant_genes/{}-{}_allp_S_{}_z_{}_b_{}{}.tsv".format(wildcards.dataset,z_col,wildcards.pinS,wildcards.pinz,wildcards.bound,suff)]
  for name in all_groups[wildcards.dataset]:
    ins.append("scripts/output/significant_genes/{}-{}_allp_S_{}_z_{}_b_{}{}.tsv".format(name,z_col,wildcards.pinS,wildcards.pinz,wildcards.bound,suff))
  print("INS",", ".join(ins))
#  return ", ".join(ins)
  return ins

def get_group_list(wildcards):
    return " ".join(all_groups[wildcards.dataset])

def get_infile(wildcards):
  return "data/{}.pq".format(wildcards.dataset)
#  if datasets[wildcards.dataset][0] == "10x":
#    return "data/{}.pq".format(wildcards.dataset)
#
#  elif datasets[wildcards.dataset][0] == "ss2":
#    return "data/{}.pq".format(wildcards.dataset)

def get_all_infiles(datasets):
  outputs = []
  for dataset in datasets:
    if datasets[dataset][0] == "10x":
      outputs.append("data/{}.pq".format(dataset))
  
    elif datasets[dataset][0] == "ss2":
      outputs.append("data/{}.pq".format(dataset))
  return outputs
print(get_SVD(datasets))

print(expand("scripts/output/variance_adjusted_permutations/{dataset}_pval-sep-ontology-sep-tiss_comp-sep-" + str(num_perms) + "_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",dataset=datasets,pinS=pins_S,pinz=pins_z,bound=bounds))
rule all:         
  input:
#    get_rijk_zscores(datasets),
#    get_SVD(datasets),
#    expand("scripts/output/perm_pvals/{dataset}_fdr_" + str(num_perms) + "_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",dataset=datasets.keys(),pinS=pins_S,pinz=pins_z,bound=bounds),
    expand("scripts/output/variance_adjusted_permutations/{dataset}_pvals_ontology-tiss_comp_" + str(num_perms) + "_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",dataset=datasets,pinS=pins_S,pinz=pins_z,bound=bounds),
    expand("scripts/output/variance_adjusted_permutations/{dataset}_pvals_compartment-tissue_" + str(num_perms) + "_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",dataset=datasets,pinS=pins_S,pinz=pins_z,bound=bounds),


#    expand("scripts/output/significant_genes/{dataset}-{z_col}_allp_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",dataset=datasets.keys(),pinS=pins_S,pinz=pins_z,bound=bounds,z_col=["scZ"])
#    get_sig(datasets, bounds),
#    get_anova(datasets),
#    get_FDR(datasets)
#
rule pq_to_tsv_SVD:
  input:
    "scripts/output/rijk_zscore/{dataset}_sym_SVD_normdonor_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".pq"

  output:
    "scripts/output/rijk_zscore/{dataset}_sym_SVD_normdonor_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv"

  resources:
    mem_mb=lambda wildcards, attempt: attempt * 40000,
#    mem_mb=lambda wildcards, attempt: attempt * 120000,

    time_min=lambda wildcards, attempt: attempt * 60
  log:
    out="job_output/pq2tsv_{dataset}_{pinS}_{pinz}_{bound}.out",
    err="job_output/pq2tsv_{dataset}_{pinS}_{pinz}_{bound}.err"

  shell:
    """
    python3.6 -u scripts/parquet_to_tsv.py --parquet {input} --outname {output} 1>> {log.out} 2>> {log.err}
    """

rule pq_to_tsv:
  input:
    "scripts/output/rijk_zscore/{dataset}_sym_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".pq"

  output:
    "scripts/output/rijk_zscore/{dataset}_sym_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv"

  resources:
    mem_mb=lambda wildcards, attempt: attempt * 40000,
#    mem_mb=lambda wildcards, attempt: attempt * 120000,

    time_min=lambda wildcards, attempt: attempt * 60
  log:
    out="job_output/pq2tsv_{dataset}_{pinS}_{pinz}_{bound}.out",
    err="job_output/pq2tsv_{dataset}_{pinS}_{pinz}_{bound}.err"

  shell:
    """
    python3.6 -u scripts/parquet_to_tsv.py --parquet {input} --outname {output} 1>> {log.out} 2>> {log.err}
    """

rule FDR_anova:
  input:
    "scripts/output/anova_zscore/{dataset}_" + z_col + "_anova_out_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",
    "scripts/output/anova_zscore/{dataset}_" + z_col + "_anova_coeff_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",


  output:
    "scripts/output/final_FDRs_anova/{dataset}_FDR_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",
    "scripts/output/final_FDRs_anova_factor/{dataset}_FDR_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv"



#    "/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/rijk_zscore/{dataset}_sym_S_{pinS}_z_{pinz}_b_{bound}.pq"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 60000,
#    mem_mb=lambda wildcards, attempt: attempt * 120000,

    time_min=lambda wildcards, attempt: attempt * 60 * 2
  log:
    out="job_output/FDR_anova_{dataset}_{pinS}_{pinz}_{bound}.out",
    err="job_output/FDR_anova_{dataset}_{pinS}_{pinz}_{bound}.err"
  params:
    suffix="_S_{pinS}_z_{pinz}_b_{bound}" + suff,
    all_datanames= get_group_list

  shell:
    """
    python3.6 -u scripts/final_FDRs_anova.py  --dataname {wildcards.dataset} --suffix {params.suffix} --all_datanames {params.all_datanames} 1>> {log.out}  2>> {log.err}
    python3.6 -u scripts/final_FDRs_anova_factor.py  --dataname {wildcards.dataset} --suffix {params.suffix} --all_datanames {params.all_datanames} 1>> {log.out}  2>> {log.err}

    """

rule significance:
  input:
    "scripts/output/rijk_zscore/{dataset}_sym_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".pq"

  output:
    "scripts/output/significant_genes/{dataset}-{z_col}_allp_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 40000,
#    mem_mb=lambda wildcards, attempt: attempt * 120000,

    time_min=lambda wildcards, attempt: attempt * 1 * 60
  log:
    out="job_output/significance_{dataset}_{pinS}_{pinz}_{bound}_{z_col}.out",
    err="job_output/significance_{dataset}_{pinS}_{pinz}_{bound}_{z_col}.err"

  params:
    unfilt={True : "--unfilt", False : ""}[unfilt],
    z_col=z_col

  shell:
    """
    python3.6 -u scripts/significant_genes.py --dataname {wildcards.dataset} --z_col {params.z_col} --pinning_S {wildcards.pinS} --pinning_z {wildcards.pinz} --lower_bound {wildcards.bound}  {params.unfilt} 1>> {log.out} 2>> {log.err}
    """

rule FDR_mz:
  input:
    "scripts/output/significant_genes/{dataset}-" + z_col + "_allp_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",
    #"/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/significant_genes/TSP2_10x_rerun_3prime-scaled_z_allp_S_0.1_z_0.0_b_5.tsv"# "/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/significant_genes/TSP1_SS2-scaled_z_allp_S_0.1_z_0.0_b_5.tsv"
#    get_group

  output:
    "scripts/output/final_FDRs_mz/{dataset}_FDR_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",

  resources:
    mem_mb=lambda wildcards, attempt: attempt * 60000,
#    mem_mb=lambda wildcards, attempt: attempt * 120000,

    time_min=lambda wildcards, attempt: attempt * 60 * 2
  log:
    out="job_output/FDR_mz_{dataset}_{pinS}_{pinz}_{bound}.out",
    err="job_output/FDR_mz_{dataset}_{pinS}_{pinz}_{bound}.err"
  params:
    suffix="_S_{pinS}_z_{pinz}_b_{bound}" + suff,
    all_datanames= get_group_list

  shell:
    """
    python3.6 -u scripts/final_FDRs_mz.py  --dataname {wildcards.dataset} --suffix {params.suffix} --all_datanames {params.all_datanames} 1>> {log.out}  2>> {log.err}
    """

rule anova:
  input:
    "scripts/output/significant_genes/{dataset}-" + z_col + "_allp_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv"

  output:
    "scripts/output/anova_zscore/{dataset}_" + z_col + "_anova_out_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",
    "scripts/output/anova_zscore/{dataset}_" + z_col + "_anova_coeff_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",
    "scripts/output/anova_factor/{dataset}_" + z_col + "_anova_out_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",
    "scripts/output/anova_factor/{dataset}_" + z_col + "_anova_coeff_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv"

  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000,
#    mem_mb=lambda wildcards, attempt: attempt * 120000,

    time_min=lambda wildcards, attempt: attempt * 5
  log:
    out="job_output/anova_{dataset}_{pinS}_{pinz}_{bound}.out",
    err="job_output/anova_{dataset}_{pinS}_{pinz}_{bound}.err"
  params:
    suffix="_S_{pinS}_z_{pinz}_b_{bound}" + suff

  shell:
    """
    Rscript scripts/anova_zscore.r  {params.suffix} {wildcards.dataset} 1>> {log.out}  2>> {log.err}
    Rscript scripts/anova_factor.r  {params.suffix} {wildcards.dataset} 1>> {log.out}  2>> {log.err}

    """


rule rijk_zscore:
  input:
    get_infile
  output:
    "scripts/output/rijk_zscore/{dataset}_sym_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".pq"
#    "/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/rijk_zscore/{dataset}_sym_S_{pinS}_z_{pinz}_b_{bound}.pq"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 80000,
#    mem_mb=lambda wildcards, attempt: attempt * 300000,

    time_min=lambda wildcards, attempt: attempt * 60 * 2 
#    time_min=lambda wildcards, attempt: attempt * 60 *  10

  log:
    out="job_output/rijk_zscore_{dataset}_{pinS}_{pinz}_{bound}.out",
    err="job_output/rijk_zscore_{dataset}_{pinS}_{pinz}_{bound}.err"

  params:
    verbose=verbose,
    light={True : "--light", False : ""}[light],
    unfilt={True : "--unfilt", False : ""}[unfilt],
    ver = ver

  shell:
    """
    python3.6 -u scripts/rijk_zscore.py {params.ver} --pinning_S {wildcards.pinS} --pinning_z {wildcards.pinz} --dataname {wildcards.dataset} --parquet {input} --lower_bound {wildcards.bound} {params.verbose} {params.light} {params.unfilt} 1>> {log.out} 2>> {log.err}
    """

rule SVD_zscore:
  input:
    "scripts/output/rijk_zscore/{dataset}_sym_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".pq"

  output:
    "scripts/output/rijk_zscore/{dataset}_sym_SVD_normdonor_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".pq"
#    "/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/rijk_zscore/{dataset}_sym_S_{pinS}_z_{pinz}_b_{bound}.pq"
  resources:
#    mem_mb=lambda wildcards, attempt: attempt * 750000,
    mem_mb=lambda wildcards, attempt: attempt * 100000,

    time_min=lambda wildcards, attempt: attempt * 60 * 6
  log:
    out="job_output/SVD_zscore_{dataset}_{pinS}_{pinz}_{bound}.out",
    err="job_output/SVD_zscore_{dataset}_{pinS}_{pinz}_{bound}.err"

  params:
    verbose=verbose,
    light={True : "--light", False : ""}[light],
    unfilt={True : "--unfilt", False : ""}[unfilt],
    ver = ver

  shell:
    """
    python3.6 -u scripts/SVD_zscore.py {params.ver} --svd_type normdonor --pinning_S {wildcards.pinS} --pinning_z {wildcards.pinz} --dataname {wildcards.dataset}  --lower_bound {wildcards.bound} {params.verbose} {params.light} {params.unfilt} 1>> {log.out} 2>> {log.err}
    """

rule perm_pval:
  input:
    "scripts/output/rijk_zscore/{dataset}_sym_SVD_normdonor_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".pq"

  output:
    "scripts/output/perm_pvals/{dataset}_fdr_" + str(num_perms) + "_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv"
#    "/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/rijk_zscore/{dataset}_sym_S_{pinS}_z_{pinz}_b_{bound}.pq"
  resources:
#    mem_mb=lambda wildcards, attempt: attempt * 750000,
    mem_mb=lambda wildcards, attempt: attempt * 60000,

    time_min=lambda wildcards, attempt: attempt * 60 * 6
  log:
    out="job_output/perm_pval_{dataset}_{pinS}_{pinz}_{bound}.out",
    err="job_output/perm_pval_{dataset}_{pinS}_{pinz}_{bound}.err"

  params:
    suffix="_S_{pinS}_z_{pinz}_b_{bound}" + suff,
    num_perms=num_perms

  shell:
    """
    python3.6 -u scripts/perm_pvals.py --suffix {params.suffix} --dataname {wildcards.dataset} --num_perms {params.num_perms} --z_col scZ 1>> {log.out} 2>> {log.err}
    python3.6 -u scripts/perm_pvals.py --z_col svd_z0 --suffix {params.suffix} --dataname {wildcards.dataset} --num_perms {params.num_perms} 1>> {log.out} 2>> {log.err}
    python3.6 -u scripts/perm_pvals.py --z_col svd_z1 --suffix {params.suffix} --dataname {wildcards.dataset} --num_perms {params.num_perms} 1>> {log.out} 2>> {log.err}
    python3.6 -u scripts/perm_pvals.py --z_col svd_z2 --suffix {params.suffix} --dataname {wildcards.dataset} --num_perms {params.num_perms} 1>> {log.out} 2>> {log.err}

    """

rule var_adj_perm_pval:
  input:
    "scripts/output/rijk_zscore/{dataset}_sym_SVD_normdonor_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".pq"

  output:
    "scripts/output/variance_adjusted_permutations/{dataset}_pvals_" + str(num_perms) + "_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",
    "scripts/output/variance_adjusted_permutations/{dataset}_outdf_" + str(num_perms) + "_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv"

  resources:
    mem_mb=lambda wildcards, attempt: attempt * 20000,

    time_min=lambda wildcards, attempt: attempt * 60 * 20
  log:
    out="job_output/var_adj_perm_pval_{dataset}_{pinS}_{pinz}_{bound}.out",
    err="job_output/var_adjperm_pval_{dataset}_{pinS}_{pinz}_{bound}.err"

  params:
    suffix="_S_{pinS}_z_{pinz}_b_{bound}" + suff,
    num_perms=num_perms

  shell:
    """
    python3.6 -u scripts/variance_adjusted_permutations.py --suffix {params.suffix} --dataname {wildcards.dataset} --num_perms {params.num_perms}  1>> {log.out} 2>> {log.err}


    """

rule var_adj_perm_pval_bytiss:
  input:
    "scripts/output/rijk_zscore/{dataset}_sym_SVD_normdonor_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".pq"

  output:
    "scripts/output/variance_adjusted_permutations/{dataset}_pvals_{group}-{sub}_" + str(num_perms) + "_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv",
    "scripts/output/variance_adjusted_permutations/{dataset}_outdf_{group}-{sub}_" + str(num_perms) + "_S_{pinS}_z_{pinz}_b_{bound}" + suff + ".tsv"

  resources:
    mem_mb=lambda wildcards, attempt: attempt * 20000,

    time_min=lambda wildcards, attempt: attempt * 60 * 20
  log:
    out="job_output/var_adj_perm_pval_bytiss_{dataset}_{pinS}_{pinz}_{bound}_{group}_{sub}.out",
    err="job_output/var_adj_perm_pval_bytiss_{dataset}_{pinS}_{pinz}_{bound}_{group}_{sub}.err"

#  wildcard_constraints:
#    group="ontology"
   
  params:
    suffix="_S_{pinS}_z_{pinz}_b_{bound}" + suff,
    num_perms=num_perms

  shell:
    """
    python3.6 -u scripts/variance_adjusted_permutations_bytiss.py --suffix {params.suffix} --dataname {wildcards.dataset} --num_perms {params.num_perms}  --group_col {wildcards.group} --sub_col {wildcards.sub} 1>> {log.out} 2>> {log.err}


    """
