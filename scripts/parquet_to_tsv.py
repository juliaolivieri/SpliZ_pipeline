import argparse
import pyarrow
import pandas as pd

def get_args():
  parser = argparse.ArgumentParser(description="convert parquet to tsv")
  parser.add_argument("-p","--parquet",help="input parquet file")
  parser.add_argument("-o","--outname",help="name to save tsv")
  parser.add_argument("--reverse",action="store_true",help="convert from tsv to pq instead")
  args = parser.parse_args()
  return args

def main():
  args = get_args()
  if args.reverse:
    df = pd.read_csv(args.outname, sep = "\t")
    df.to_parquet(args.parquet)
  else:
    df = pd.read_parquet(args.parquet)
    df.to_csv(args.outname, sep = "\t", index = False)
  
main()
