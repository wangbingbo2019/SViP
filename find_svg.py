import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from pandas.core.frame import DataFrame
import scanpy as sc
import scipy as sp
import matplotlib.image as mpimg
from matplotlib import rcParams
import seaborn as sb
import anndata as ad
import h5py
from PIL import Image
import cv2
import sys
import SpatialDE
import NaiveDE
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

'''
count_csv: spot_id * gene_id
location_csv: need columns: x, y
'''
def find_svg(count_csv, location_csv, out_excel, method='SpatialDE'):
    if method == 'SpatialDE':
        df_count = pd.read_csv(count_csv, sep='\t',index_col=0)
        counts = df_count.T[df_count.sum(0) >= 3].T

        loc=pd.read_csv(location_csv, sep='\t',index_col=0)
        loc1=loc.loc[counts.index,'x'].values
        loc2=loc.loc[counts.index,'y'].values

        sample_info = pd.DataFrame(columns = ['','x','y','total_counts'])


        sample_info['']=counts.index.tolist()
        sample_info['x']=loc1
        sample_info['y']=loc2
        sample_info['total_counts']=counts.T.sum(0).tolist()
        sample_info = sample_info.set_index('')

        norm_expr = NaiveDE.stabilize(counts.T).T
        norm_expr = counts

        del counts
        # norm_expr = NaiveDE.regress_out(sample_info,norm_expr.T,'np.log(total_counts)').T
        print(norm_expr)

        X = sample_info[['x','y']]

        results = SpatialDE.run(X,norm_expr)
        res = results.sort_values('qval')[['g','pval','qval']]
        res.to_excel(out_excel, index=False)
    elif method == 'SPARK':
        print('using SPARK\'s SVG...')
    else:
        print('Please give an avaliable method:[SpatialDE, SPARK]')












