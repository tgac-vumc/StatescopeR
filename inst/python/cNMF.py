#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BLADE.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Run cNMF for state discovery
# *Edited from Jurriaan Janssens (j.janssen4@amsterdamumc.nl) version for use 
# in R
#
# Author: Mischa Steketee (m.f.b.steketee@amsterdamumc.nl)
#
#
# TODO:
# 1) 
#
# History:
#  09-12-2024: Edited for R from Python version
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
import sys
import cNMF_functions
import numpy as np
import pandas as pd
import random
from scipy.stats import zscore
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_score
import random
from scipy.cluster.hierarchy import average, cophenet
from scipy.spatial.distance import  pdist
from joblib import Parallel, delayed
import dist

#-------------------------------------------------------------------------------
# 1.1 Run initial NMF runs
#------------------------------------------------------------------------------- 
data_dict = dict()
for k in range(2,args.max_clusters):
    print(k)
    cNMF_model, cophcor, consensus_matrix = cNMF(data_scaled, k, args.n_iter, args.Ncores)
    H = cNMF_model.H
    cluster_assignments = []
    for i in range(H.shape[1]):
        cluster_assignments.append(int(np.where(H[:,i] == max(H[:,i]))[0] + 1))    
    data_dict[k] = {'model':cNMF_model,'cophcor':cophcor, 'consensus': consensus_matrix,'cluster_assignments':cluster_assignments}

#-------------------------------------------------------------------------------
# 3.2 Choose k
#------------------------------------------------------------------------------- 
cophcors = [d['cophcor'] for d in data_dict.values()]
ks = [k for k in data_dict.keys()]

nclust = find_threshold(cophcors,ks,min_cophenetic=args.min_coph)
drop = biggest_drop(cophcors)
if not nclust:
    nclust = drop


# Plot cophenetic coefficients
fig, ax = plt.subplots(figsize = (5,3))
ax.plot(ks,cophcors, c = 'black',zorder =-1)
for cophcor,k in zip(cophcors,ks):
    if k == nclust:
        ax.scatter(k,cophcor,c = 'red')
    else:
        ax.scatter(k,cophcor, c= 'black')

ax.set_xticks(ks,ks)
ax.set_xlabel('k')
plt.tight_layout()
fig.savefig(args.output_coph)

#-------------------------------------------------------------------------------
# 3.2 Run final model
#------------------------------------------------------------------------------- 
print('Running final model:')
# Run final model    
cNMF_model, cophcor, consensus_matrix = cNMF(data_scaled, nclust, args.n_iter_final, args.Ncores)
H = cNMF_model.H





cluster_assignments = []
for i in range(H.shape[1]):
    cluster_assignments.append(int(np.where(H[:,i] == max(H[:,i]))[0] + 1))    
data_dict['final'] = {'model':cNMF_model,'cophcor':cophcor, 'consensus': consensus_matrix,'cluster_assignments':cluster_assignments}

# Store input matrix
data_dict['input_matrix'] = data_scaled

# Fetch output files
H = pd.DataFrame(cNMF_model.H.transpose(), index = Purification_matrix.index)
W = pd.DataFrame(cNMF_model.W, index = Purification_matrix.columns)
