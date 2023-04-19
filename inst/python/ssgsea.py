# import
import scanpy as sc
import pandas as pd
import gseapy as gp
import argparse


# Parse the parameters passed in
parser = argparse.ArgumentParser(description=' ')
parser.add_argument('--min_size', type=int, default = None)
parser.add_argument('--max_size', type=int, default = None)
parser.add_argument('--seed', type=int, default = None)
parser.add_argument('--threads', type=int, default = None)

args = parser.parse_args()

# assign
min_size = args.min_size
max_size = args.max_size
seed = args.seed
threads = args.threads

# def function
adata = sc.read_h5ad("./temp.h5ad")
data = pd.DataFrame.sparse.from_spmatrix(adata.X)
data.index=adata.obs.index
data.columns=adata.var.index
data = data.T

geneset = pd.read_csv("./geneset.csv")

gene_sets = geneset.groupby('source')['target'].apply(list).to_dict()

ss = gp.ssgsea(data= data,
               gene_sets=gene_sets,
               outdir=None,
               sample_norm_method='rank', # choose 'custom' will only use the raw value of `data`
               no_plot=True,
               threads=threads,
               seed=seed,
               min_size=min_size,
               max_size=max_size)

nes = ss.res2d.pivot(index='Term', columns='Name', values='ES')
# output result
nes.to_csv('./matrix.py.result.csv')








