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
import numpy as np
import pandas as pd
import random
from scipy.stats import zscore
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import random
from scipy.cluster.hierarchy import average, cophenet
from scipy.spatial.distance import  pdist
from joblib import Parallel, delayed
import logging

#----------------------------------------------------------------------------
# 1.2 Define functions for using cNMF
#------------------------------------------------------------------------------- 
## From https://github.com/jurriaanjanssen/scRNAseq-snakemake
# Function to extract hard cluster assignments form soft clusters
def cluster_assignment(H):
    return [np.argmax(H[:,i]) for i in range(H.shape[1])]

# Function to calculate connectivity matrix for hard cluster assignments
def connectivity(clusters):
    Nsample = len(clusters)
    connectivity_matrix = np.zeros((Nsample,Nsample))
    for i in range(Nsample):
        connectivity_matrix[i,np.where(clusters == clusters[i])] = 1
    return connectivity_matrix

# calculate cophenetic coefficient for consensus matrices
def Calculate_cophcor(X):
    dmat = pdist(X)
    Z = average(dmat)
    coph = cophenet(Z)
    return np.corrcoef(dmat,coph)[0,1]

    
def cNMF(data,k,nrun,ncores,niter=1000):
    def Run_cNMF(data,k,seed,niter=1000):
        random.seed(seed)
        model = CNMF(data, num_bases=k, niter=niter)
        model.factorize(niter=niter)
        return model
    
    # perform NMF in parralel with different random seeds
    models = Parallel(n_jobs=ncores, verbose=10)(
        delayed(Run_cNMF)(data, k, i, niter)
        for i in range(nrun))
    # Fetch cluster assignments, connectivity matrices and calculate the consensus matrix
    consensus_matrix = sum([connectivity(cluster_assignment(mod.H)) for mod in models]) / nrun
    # Calculate the cophenetic correlation coefficient
    cophcor = Calculate_cophcor(consensus_matrix)
    # Calculate objective function and select best model
    objective = [mod.frobenius_norm() for mod in models]
    # pick model with lowest objective value
    model = models[objective.index(min(objective))]
    return model, cophcor, consensus_matrix

                
def find_threshold(cophcors,ks,min_cophenetic=0.95):
    cross = [False, False] + [((cophcors[i] < min_cophenetic) & (cophcors[i-1] >= min_cophenetic) & (cophcors[i-2] >= min_cophenetic)) for i in range(2,len(cophcors))]
    if any(cross) == False:
        return ks[-1]
    else:
        last_cross = [i for i,x in enumerate(cross) if x][-1]
        diff = [abs(min_cophenetic-cophcors[last_cross-1]),min_cophenetic-cophcors[last_cross]]
        return ks[last_cross - 2 + [i for i,x in enumerate(diff) if x == min(diff)][0]]


def biggest_drop(cophcors):
    cross = [-100] + [cophcors[i-1] - cophcors[i] for i in range(1,len(cophcors))]
    return cross.index(max(cross))-1


def select_k(data_scaled, max_clusters, n_iter, Ncores, min_cophenetic):
    data_dict = dict()
    for k in range(2, max_clusters):
        cNMF_model_k, cophcor_k, consensus_k = cNMF(data_scaled, k,  n_iter, Ncores)
        
        H = cNMF_model_k.H
        cluster_assignments = [
                int(np.where(H[:, i] == H[:, i].max())[0] + 1)
                for i in range(H.shape[1])
        ]
        data_dict[k] = {
                "model":              cNMF_model_k,
                "cophcor":            cophcor_k,
                "consensus":          consensus_k,
                "cluster_assignments": cluster_assignments,
            }
    
    ks         = sorted(data_dict)
    sweep_curve = [data_dict[k]["cophcor"] for k in ks]
    nclust      = find_threshold(sweep_curve, ks, min_cophenetic) \
                      or biggest_drop(sweep_curve)
        
    return nclust    
                
