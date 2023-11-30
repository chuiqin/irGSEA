# import
import os
import pandas as pd 
import numpy as np
import scanpy as sc
import decoupler as dc
import argparse
from scipy.sparse import csr_matrix

# Parse the parameters passed in
parser = argparse.ArgumentParser(description=' ')
parser.add_argument('--kcdf', type=str, default = None)
parser.add_argument('--min_n', type=int, default = None)
parser.add_argument('--seed', type=int, default = None)

args = parser.parse_args()

# assign
kcdf = args.kcdf
min_n = args.min_n
seed = args.seed



# def function
adata = sc.read_h5ad("./temp.h5ad")
geneset = pd.read_csv("./geneset.csv")

dc.run_viper(
    mat=adata,
    net=geneset,
    source='source',
    target='target',
    verbose=True,
    use_raw=False
)


# output result
adata.obsm['viper_estimate'].to_csv('./matrix.py.result.csv')

