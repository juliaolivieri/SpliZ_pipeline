import argparse
import datetime 
import numpy as np
import pandas as pd
import pickle
import pyarrow
import time
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

def get_args():
  parser = argparse.ArgumentParser(description="calculate splicing scores per gene/cell")
  parser.add_argument("--dataname",help="name of dataset to use")
  parser.add_argument("--parquet",help="input parquet file")
  parser.add_argument("--pinning_S",type=float,help="pinning level for S_ijks")
  parser.add_argument("--pinning_z",type=float,help="pinning level for zs")
  parser.add_argument("--light",action="store_true",help="if included, don't calculate extra columns (saves time)")
  parser.add_argument("--verbose",action="store_true",help="print times")
  parser.add_argument("--unfilt",action="store_true",help="don't filter by SICILIAN")
  parser.add_argument("--v2",action="store_true",help="filter by SICILIAN v2")


#  parser.add_argument("--bound_lower",action="store_true",help="include lower bound on number of junctional reads a cell + gene needs to have in order to get a z score")

  parser.add_argument("--lower_bound",type=int,help="only include cell/gene pairs the have at least this many junctional reads for the gene")


  args = parser.parse_args()
  return args

def prepare_df(df, let, rank_by_donor, rev_let, let_dict):
  
  # create donor identifier
  df["pos{}_group".format(let)] = df["juncPosR1{}".format(let)].astype(str) + df["geneR1A_uniq"]
  df["rank_" + let_dict[let]] = df.groupby("pos{}_group".format(let))["juncPosR1{}".format(rev_let[let])].rank(method="dense")

  # remove consitutive splicing
  df["max_rank"] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["rank_" + let_dict[let]].max())
  df = df[df["max_rank"] > 1]
  
  if not rank_by_donor:
    df["rank_" + let_dict[let]] = df.groupby("geneR1A_uniq")["juncPosR1B"].rank(method="dense")
  return df

def calc_Sijk(df,let, pinning_S, let_dict):
  # calculate the average rank calculation per gene
  # same as this calculation (this one's slower): df["rank_mean"] = df.groupby("pos{}_group".format(let)).apply(lambda x: (x["numReads"] * x["rank_acc"])/x["numReads"].sum()).reset_index(level=0,drop=True)

  # number of reads with this donor across all cells
  df["sum_reads_group"] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["numReads"].sum())

  df["read_x_" + let_dict[let]] = df["numReads"] * df["rank_" + let_dict[let]]

  # the sum of acceptors for all reads in all cells with this donor
  df["num"] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["read_x_" + let_dict[let]].sum())

  # average acceptor for a read with this donor (donor has one value for this)
  df["rank_mean"]= df["num"] / df["sum_reads_group"]

#   df["rank_mean"] = df.groupby("pos{}_group".format(let)).apply(lambda x: (x["numReads"] * x["rank_acc"])/x["numReads"].sum())
  # sum squared difference in rank for ever read
  df["sq_diff"] = df["numReads"] * (df["rank_" + let_dict[let]] - df["rank_mean"])**2

  # Get the sum of these squared differences for each donor
  df["don_num"] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["sq_diff"].sum())

  # sum of squared differences normalized by total number of reads
  # changed to make it the sample standard deviation (added minus 1)
  df["don_sigma"] = df["don_num"] / (df["sum_reads_group"])

#   df["don_sigma"] = df.groupby("pos_{}_group".format(let)).apply(lambda x: (x["numReads"] * (x["rank_acc"] - x["rank_mean"])**2).sum()/x["numReads"].sum())
  # this is the S_ijk value (difference normalized by sd) - should be normal 0/1
  df["S_ijk_{}".format(let)] = (df["rank_" + let_dict[let]] - df["rank_mean"])/np.sqrt(df["don_sigma"])

  # round outlying S values
  low_quant = df["S_ijk_{}".format(let)].quantile(q=pinning_S)
  high_quant = df["S_ijk_{}".format(let)].quantile(q=1 - pinning_S)
  df["S_ijk_{}_unpinned".format(let)] = df["S_ijk_{}".format(let)]

  df.loc[df["S_ijk_{}".format(let)] < low_quant,"S_ijk_{}".format(let)] = low_quant
  df.loc[df["S_ijk_{}".format(let)] > high_quant,"S_ijk_{}".format(let)] = high_quant



  # correct for those with no variance
  df.loc[df["don_sigma"] == 0, "S_ijk_{}".format(let)] = 0
  df["n_sijk"] = df["numReads"]
  df.loc[df["don_sigma"] == 0,"n_sijk"] = 0
  return df

def normalize_Sijks(df,let):
  # calculate mean of SijkA's per gene
  df["n_s"] = df["numReads"] * df["S_ijk_" + let]
  df["num"] = df["geneR1A_uniq"].map(df.groupby("geneR1A_uniq")["n_s"].sum())
  df["n_gene"] = df["geneR1A_uniq"].map(df.groupby("geneR1A_uniq")["numReads"].sum())
  df["sijk{}_mean".format(let)] = df["num"] / df["n_gene"]

  # calculate standard deviation of SijkA's per gene
  # temp_df["n_g"] = temp_df["cell_gene"].map(temp_df.groupby("cell_gene")["n_sijk"].sum())
  # temp_df["sijkA_mean"] = temp_df["num"] / temp_df["n_gene"]
  # temp_df["S_ijk_A_mean"] = temp_df["geneR1A_uniq"].map(temp_df.groupby("geneR1A_uniq")["S_ijk_A"].mean())
  df["sd_num"] = df["numReads"] * (df["S_ijk_" + let] - df["sijk{}_mean".format(let)])**2
  df["num"] = df["geneR1A_uniq"].map(df.groupby("geneR1A_uniq")["sd_num"].sum())
  df["sijk{}_var".format(let)] = df["num"] / df["n_gene"]
#  df["S_ijk_{}".format(let)] = (df["S_ijk_{}".format(let)] - df["sijk{}_mean".format(let)])/np.sqrt(df["sijk{}_var".format(let))
  return df

def main():
  SVD = False
  let_dict = {"A" : "acc", "B" : "don"}
  t0 = time.time()
  args = get_args()
  outpath = "scripts/output/rijk_zscore/"

  df = pd.read_parquet(args.parquet,columns=["inc_emp.p","tissue","gene_count_per_cell_filt","juncPosR1A","geneR1A_uniq","juncPosR1B","numReads","cell","channel","splice_ann","tissue","compartment","free_annotation","missing_domains","domain_insertions","refName_newR1","gene_frac_filt","called","chrR1A","exon_annR1A","exon_annR1B"])
  if args.verbose:
    print("read in parquet",datetime.timedelta(seconds = time.time() - t0))
  if "missing_domains" in df.columns and not args.light:
    domain_breakdown = True
  else:
    domain_breakdown = False

  df.reset_index(drop=True,inplace=True)
  rank_by_donor = True
#  no_panc = True

#  bound_lower = True
#  lower_gene_count = 10
#  pinning = 0.0001


  if args.unfilt:

    # only include junctions with more than 1 read in the dataset
    df["numReads_tot"] = df["refName_newR1"].map(df.groupby("refName_newR1")["numReads"].sum())
    df = df[df["numReads_tot"] > 1]

  elif args.v2:
    print("SICILIAN v2")
    df = df[df["called"] == 1]

  else:

    # only SICILIAN junctions
    df = df[df["inc_emp.p"]]

  
#  df["cell_gene"] = df["cell"] + df["geneR1A_uniq"]


 # print("RPL29 shape B",df[df["geneR1A_uniq"] == "RPL29"].shape)

  # use second location gene name if first is unknown
  df["geneR1B_uniq"] = df["refName_newR1"].str.split("|").str[1].str.split(":").str[1]
  idx = df[(df["geneR1A_uniq"].isin(["unknown",""])) | (df["geneR1A_uniq"].isna())].index
  df.loc[idx,"geneR1A_uniq"] = df.loc[idx,"geneR1B_uniq"]

  bin_size = 100000
  # bin unknown genes
  idx = df[(df["geneR1A_uniq"] == "") | (df["geneR1A_uniq"] == "unknown") | (df["geneR1A_uniq"].isna())].index
  df.loc[idx,"geneR1A_uniq"] = "unknown_" + df["chrR1A"].astype(str) + "_" + (df.loc[idx]["juncPosR1A"] - df.loc[idx]["juncPosR1A"] % bin_size).astype(str)
  print("replaced gene names",df[(df["geneR1A_uniq"].isin(["unknown",""])) | (df["geneR1A_uniq"].isna())].shape[0])

  if args.verbose:
    print("replace with geneR1B",datetime.timedelta(seconds = time.time() - t0))

  # get sign of gene to adjust z score
  sign_df = df.drop_duplicates("geneR1A_uniq")
  sign_df["strandA"] = sign_df["refName_newR1"].str.split("|").str[0].str.split(":").str[3]
  sign_df["strandB"] = sign_df["refName_newR1"].str.split("|").str[1].str.split(":").str[3]
  idx = sign_df[sign_df["strandA"] == "?"].index
  sign_df.loc[idx,"strandA"] = sign_df.loc[idx,"strandB"]
  sign_df["sign"] = 1
  sign_df.loc[sign_df["strandA"] == "-","sign"] = -1
  sign_df[["geneR1A_uniq","strandA","sign"]]
  sign_dict = pd.Series(sign_df.sign.values,index=sign_df.geneR1A_uniq).to_dict()
  df["sign"] = df["geneR1A_uniq"].map(sign_dict).fillna(1)
  if args.verbose:
    print("get sign",datetime.timedelta(seconds = time.time() - t0))

#  df = df[df["tissue"] != "Pancreas"]

  # only cells we have annotations for
#  df = df[~df["free_annotation"].isna()]

  # remove low quality cells
#  df["cell_counts"] = df["cell"].map(df["cell"].value_counts())
#  df = df[df["cell_counts"] > 200]
  df["cell_gene"] = df["cell"] + df["geneR1A_uniq"]

  rev_let = {"A" : "B", "B" : "A"}

  if domain_breakdown:
    split_dict = {True : ["ann", "dom_ch"], False : ["unann", "dom_unch"]}
  else:
    split_dict = {True : ["ann"], False : ["unann"]}

  # remove constitutive splicing
  df["posA_group"] = df["juncPosR1A"].astype(str) + df["geneR1A_uniq"]
  df["posB_group"] = df["juncPosR1B"].astype(str) + df["geneR1A_uniq"]

  df["rank_acc"] = df.groupby("posA_group")["juncPosR1B"].rank(method="dense")
  df["rank_don"] = df.groupby("posB_group")["juncPosR1A"].rank(method="dense")

  df["max_rank_acc"] = df["posA_group"].map(df.groupby("posA_group")["rank_acc"].max())
  df["max_rank_don"] = df["posB_group"].map(df.groupby("posB_group")["rank_don"].max())

  # add domain columns
  letters = ["A","B"]
  for let in letters:

    if domain_breakdown:
      df["num_missing_" + let] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["missing_domains"].nunique())
      df["num_inserted_" + let] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["domain_insertions"].nunique())
      df["domain_changed_" + let] = (df["num_missing_" + let] + df["num_inserted_" + let]) > 0


  df = df[(df["max_rank_don"] > 1) | (df["max_rank_acc"] > 1)]
  if args.verbose:
    print("remove constitutive",datetime.timedelta(seconds = time.time() - t0))
  # require at least args.lower_bound nonconstitutive spliced reads
  df["noncon_count"] = df["cell_gene"].map(df.groupby("cell_gene")["numReads"].sum())
  df = df[df["noncon_count"] > args.lower_bound]

  full_df = df.copy()

  calc_dfs = {}

  for let in tqdm(["A","B"]):
    df = full_df
    # create donor identifier
    df = prepare_df(df, let, rank_by_donor, rev_let, let_dict)

    if args.verbose:
      print("prepare_df",datetime.timedelta(seconds = time.time() - t0))
    df = calc_Sijk(df,let,args.pinning_S, let_dict)
    if args.verbose:
      print("calc_sijk",datetime.timedelta(seconds = time.time() - t0))

    df = normalize_Sijks(df,let)
    if args.verbose:
      print("normalize_Sijk",datetime.timedelta(seconds = time.time() - t0))

    # remove those with variance == 0
    df = df[df["sijk{}_var".format(let)] != 0]

    # calculate z score 
#    df["n_g"] = df["cell_gene"].map(df.groupby("cell_gene")["n_sijk"].sum())
    df["n.g_" + let] = df["cell_gene"].map(df.groupby("cell_gene")["numReads"].sum())

    df["nSijk" + let] = (df["S_ijk_" + let] - df["sijk{}_mean".format(let)]) / np.sqrt(df["sijk{}_var".format(let)])
    df["mult"] = df["numReads"] * df["nSijk" + let]  / np.sqrt(df["n.g_" + let])
    df["z_" + let] = df["sign"] * df["cell_gene"].map(df.groupby("cell_gene")["mult"].sum())
    df["scaled_z_" + let] = df["z_" + let] / np.sqrt(df["n.g_" + let])

    if args.verbose:
      print("calc_z",datetime.timedelta(seconds = time.time() - t0))

    ############## end modify Sijk ####################
    df["cell_gene_junc"] = df["cell_gene"] + df["refName_newR1"]

    if not args.light:
      # calculate the z score
      df["x_sijk"] = df["S_ijk_{}".format(let)] * df["n_sijk"]
#      df["x_sijk"] = df["S_ijk_{}".format(let)] * df["n_sijk"]

      df["num"] = df["cell_gene"].map(df.groupby("cell_gene")["x_sijk"].sum())
      df["denom_sq"] = df["cell_gene"].map(df.groupby("cell_gene")["n_sijk"].sum())
  #    df["z_{}".format(let)] = df["sign"] * df["num"]/np.sqrt(df["denom_sq"])
  
      # get junction that "contributes the most" to the z score
      df["temp"] = df["x_sijk"] / np.sqrt(df["denom_sq"])
      df["temp_mag"] = abs(df["temp"])
      df["idxmax_z"] = df["cell_gene"].map(df.groupby("cell_gene")["temp_mag"].idxmax())
      map_df = df.loc[df["idxmax_z"],["cell_gene","refName_newR1","temp"]]
      df["junc_max_{}".format(let)] = df["cell_gene"].map(pd.Series(map_df.refName_newR1.values,index=map_df.cell_gene).to_dict())
      df["max_don_z_{}".format(let)] = df["cell_gene"].map(pd.Series(map_df.temp.values,index=map_df.cell_gene).to_dict())

#    # break down z score by annotation
#    df["num_ann"] = df["cell_gene"].map(df[df["splice_ann"]].groupby("cell_gene")["x_sijk"].sum())
#    df["z_{}_ann".format(let)] = df["sign"] * df["num_ann"]/np.sqrt(df["denom_sq"])
#
#    df["num_unann"] = df["cell_gene"].map(df[~df["splice_ann"]].groupby("cell_gene")["x_sijk"].sum())
#    df["z_{}_unann".format(let)] = df["sign"] * df["num_unann"]/np.sqrt(df["denom_sq"])
  
    if args.pinning_z != 0:
      # round outlying z values
      low_quant = df["z_{}".format(let)].quantile(q=args.pinning_z)
      high_quant = df["z_{}".format(let)].quantile(q=1 - args.pinning_z)
      
      df.loc[df["z_{}".format(let)] < low_quant,"z_{}".format(let)] = low_quant
      df.loc[df["z_{}".format(let)] > high_quant,"z_{}".format(let)] = high_quant

    if not args.light:
      # break down z score by annotation
      for k,v in split_dict.items():
  #      print(df[df["splice_ann"] == k].shape)
        df["num_{}".format(v[0])] = df["cell_gene"].map(df[df["splice_ann"] == k].groupby("cell_gene")["x_sijk"].sum())
  
        if domain_breakdown:
          df["num_{}".format(v[1])] = df["cell_gene"].map(df[df["domain_changed_" + let] == k].groupby("cell_gene")["x_sijk"].sum())
  
        for y in v:
  
          df["z_{}_{}".format(let,y)] = df["sign"] * df["num_{}".format(y)]/np.sqrt(df["denom_sq"])
      
          # round outlying z values
          low_quant = df["z_{}_{}".format(let,y)].quantile(q=args.pinning_z)
          high_quant = df["z_{}_{}".format(let,y)].quantile(q=1 - args.pinning_z)
      
          df.loc[df["z_{}_{}".format(let,y)] < low_quant,"z_{}_{}".format(let,y)] = low_quant
          df.loc[df["z_{}_{}".format(let,y)] > high_quant,"z_{}_{}".format(let,y)] = high_quant
#    df.to_parquet("{}{}_sym_S_{}_z_{}_b_{}_light_{}.pq".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound,let))
       

    calc_dfs[let] = df

  df = calc_dfs["A"].merge(calc_dfs["B"],on="cell_gene_junc",how="outer",suffixes=("","_x"))
  if args.verbose:
    print("merged")
  for cx in [x for x in df.columns if x.endswith("_x")]:
    c = cx[:-2]
    df.loc[df[c].isna(),c] = df.loc[df[c].isna(),cx]
  df.drop([x for x in df.columns if x.endswith("_x")],inplace=True,axis=1)


  # average two scores (negate one of them)

  grouped = df.groupby('geneR1A_uniq')
  for let in ["A","B"]:
    z_dict = pd.Series(calc_dfs[let]["z_" + let].values,index=calc_dfs[let].cell_gene).to_dict()
    df["z_" + let] = df["cell_gene"].map(z_dict)
    scz_dict = pd.Series(calc_dfs[let]["scaled_z_" + let].values,index=calc_dfs[let].cell_gene).to_dict()
    df["scaled_z_" + let] = df["cell_gene"].map(z_dict)

  df["cov"] = df["geneR1A_uniq"].map(grouped.apply(lambda x: x['z_A'].cov(x['z_B'])))

  idx = df[df["z_A"].isna()].index
  df.loc[idx,"z"] = -df.loc[idx,"z_B"]
  df.loc[idx,"scZ"] = -df.loc[idx,"scaled_z_B"]

  idx = df[df["z_B"].isna()].index
  df.loc[idx,"z"] = df.loc[idx,"z_A"]
  df.loc[idx,"scZ"] = df.loc[idx,"scaled_z_A"]


  idx = df[(~df["z_A"].isna()) & (~df["z_B"].isna())].index
  df.loc[idx,"z"] = (df.loc[idx,"z_A"] - df.loc[idx,"z_B"])/np.sqrt(2 )
  df.loc[idx,"scZ"] = (df.loc[idx,"scaled_z_A"] - df.loc[idx,"scaled_z_B"])/np.sqrt(2 )

  if args.verbose:
    print("avg z",datetime.timedelta(seconds = time.time() - t0))

#  df["negz_B"] = -1*df["z_B"]
#  df["z"] = df[["z_A","negz_B"]].mean(axis=1)

  #df["negz_B"] = -1*df["z_B_ann"]
  #df["z_ann"] = df[["z_A_ann","negz_B"]].mean(axis=1)

  #df["negz_B"] = -1*df["z_B_unann"]
  #df["z_unann"] = df[["z_A_unann","negz_B"]].mean(axis=1)

  if not args.light:
    # average two scores for split z
    for v in split_dict.values():
      for y in v:
        grouped = df.groupby('geneR1A_uniq')
        df["cov_{}".format(y)] = df["geneR1A_uniq"].map(grouped.apply(lambda x: x['z_A_{}'.format(y)].cov(x['z_B_{}'.format(y)])))
  
        idx = df[df["z_A_{}".format(y)].isna()].index
        df.loc[idx,"z_{}".format(y)] = -df.loc[idx,"z_B_{}".format(y)]
      
        idx = df[df["z_B_{}".format(y)].isna()].index
        df.loc[idx,"z_{}".format(y)] = df.loc[idx,"z_A_{}".format(y)]
      
        idx = df[(~df["z_A_{}".format(y)].isna()) & (~df["z_B_{}".format(y)].isna())].index
        df.loc[idx,"z_{}".format(y)] = (df.loc[idx,"z_A_{}".format(y)] - df.loc[idx,"z_B_{}".format(y)])/np.sqrt(2) - df["cov_{}".format(y)]

#      df["negz_B_{}".format(y)] = -1*df["z_B_{}".format(y)]
#      df["z_{}".format(y)] = df[["z_A_{}".format(y),"negz_B_{}".format(y)]].mean(axis=1)

#  print("RPL29 shape F",df[df["geneR1A_uniq"] == "RPL29"].shape)

#  df.drop(["denom_sq","num_ann","num_unann","num","temp","temp_mag"],axis=1,inplace=True)
#  df.drop(["negz_B_unann","negz_B_ann","negz_B","denom_sq","num_ann","num_unann","num","temp","temp_mag"],axis=1,inplace=True)

  df["ontology"] = df["tissue"] + df["compartment"] + df["free_annotation"]
  
#  df["ontology_gene"] = df["ontology"] + df["geneR1A_uniq"]
  suff = ""
  if args.light:
    suff += "_light"
  if args.unfilt:
    suff += "_unfilt"

  df["n.g"] = df["cell_gene"].map(df.groupby("cell_gene")["numReads"].sum())
  df["scaled_z"] = df["z"] / np.sqrt(df["n.g"])

  for let in ["A","B"]:
    df["zcontrib" + let] = df["numReads"] * df["nSijk" + let] / np.sqrt(df["n.g"])

  ##### PERFORM SVD ZSCORE CALCULATION #####

  if SVD:
    letters = ["A","B"]
    for let in letters:
      df["zcontrib{}_rep".format(let)] = df["zcontrib" + let].fillna(0)
      df["juncPosR1" + let] = df["juncPosR1" + let].astype(int).astype(str) + "_" + let
  
    k = 3 # number of components to include
    loads = {"f{}".format(i) : {} for i in range(k)}
    zs = {"svd_z{}".format(i) : {} for i in range(k)}
    
    for gene, gene_df in tqdm(df.groupby("geneR1A_uniq")):
      
      # get zcontrib matrix
      gene_mats = []
      for let in letters:
        gene_mat = gene_df.pivot_table(index="cell_gene",columns="juncPosR1{}".format(let),values="zcontrib{}_rep".format(let),fill_value=0)
        gene_mats.append(gene_mat)
      gene_mat = gene_mats[0].merge(gene_mats[1],on="cell_gene")
      
      # calculate svd
      u, s, vh = np.linalg.svd(gene_mat)
      
      if len(s) >= k:
        # calculate new z scores based on svd
        new_zs = gene_mat.dot(np.transpose(vh[:k,:]))
    
        # calculate load on each component
        load = np.square(s)/sum(np.square(s))
    
        # save new zs and fs in dictionaries to save later
        for i in range(k):
          loads["f{}".format(i)][gene] = load[i]
          zs["svd_z{}".format(i)].update(pd.Series(new_zs[i].values,index=new_zs.index).to_dict())
    
        # save loadings
        v_out = pd.DataFrame(vh,columns=gene_mat.columns)
        v_out.to_csv("{}SVD/{}_{}_S_{}_z_{}_b_{}{}.tsv".format(outpath, gene,args.dataname, args.pinning_S, args.pinning_z, args.lower_bound, suff), index=False, sep = "\t")
        
    for i in range(k):
      df["f{}".format(i)] = df["geneR1A_uniq"].map(loads["f{}".format(i)])
      df["svd_z{}".format(i)] = df["cell_gene"].map(zs["svd_z{}".format(i)])
    df["svd_z_sumsq"] = (df[["svd_z{}".format(i) for i in range(k)]]**2).sum(axis=1)

    sub_cols = ["cell","geneR1A_uniq","tissue","compartment","free_annotation","ontology","scZ","svd_z_sumsq","n.g_A","n.g_B"] + ["f{}".format(i) for i in range(k)] + ["svd_z{}".format(i) for i in range(k)]
  else:
    sub_cols = ["cell","geneR1A_uniq","tissue","compartment","free_annotation","ontology","scZ","n.g_A","n.g_B"] 

  df.drop_duplicates("cell_gene")[sub_cols].to_csv("{}{}_sym_S_{}_z_{}_b_{}{}_subcol.tsv".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff),index=False,sep="\t")
  df.to_parquet("{}{}_sym_S_{}_z_{}_b_{}{}.pq".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff))

  if args.verbose:
    print("wrote parquet",datetime.timedelta(seconds = time.time() - t0))

  if args.verbose:
    print("wrote files",datetime.timedelta(seconds = time.time() - t0))



main()
