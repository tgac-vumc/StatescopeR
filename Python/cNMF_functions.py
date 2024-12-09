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
from OncoBLADE import Framework_Iterative
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
import logging
import dist

#-------------------------------------------------------------------------------
# 1.1 Define functions for actual cNMF
#------------------------------------------------------------------------------- 
## From: https://github.com/cthurau/pymf

class PyMFBase():
    """
    PyMF Base Class. Does nothing useful apart from poviding some basic methods.
    """
    # some small value
   
    _EPS = np.finfo(float).eps
    
    def __init__(self, data, num_bases=4, **kwargs):
        """
        """
        
        def setup_logging():
            # create logger       
            self._logger = logging.getLogger("pymf")
       
            # add ch to logger
            if len(self._logger.handlers) < 1:
                # create console handler and set level to debug
                ch = logging.StreamHandler()
                ch.setLevel(logging.DEBUG)
                # create formatter
                formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        
                # add formatter to ch
                ch.setFormatter(formatter)

                self._logger.addHandler(ch)

        setup_logging()
        
        # set variables
        self.data = data       
        self._num_bases = num_bases             
      
        # initialize H and W to random values
        self._data_dimension, self._num_samples = self.data.shape
        

    def residual(self):
        """ Returns the residual in % of the total amount of data

        Returns
        -------
        residual : float
        """
        res = np.sum(np.abs(self.data - np.dot(self.W, self.H)))
        total = 100.0*res/np.sum(np.abs(self.data))
        return total
        
    def frobenius_norm(self):
        """ Frobenius norm (||data - WH||) of a data matrix and a low rank
        approximation given by WH. Minimizing the Fnorm ist the most common
        optimization criterion for matrix factorization methods.

        Returns:
        -------
        frobenius norm: F = ||data - WH||

        """
        # check if W and H exist
        if hasattr(self,'H') and hasattr(self,'W'):
            if scipy.sparse.issparse(self.data):
                tmp = self.data[:,:] - (self.W * self.H)
                tmp = tmp.multiply(tmp).sum()
                err = np.sqrt(tmp)
            else:
                err = np.sqrt( np.sum((self.data[:,:] - np.dot(self.W, self.H))**2 ))            
        else:
            err = None

        return err
        
    def _init_w(self):
        """ Initalize W to random values [0,1].
        """
        # add a small value, otherwise nmf and related methods get into trouble as 
        # they have difficulties recovering from zero.
        self.W = np.random.random((self._data_dimension, self._num_bases)) + 10**-4
        
    def _init_h(self):
        """ Initalize H to random values [0,1].
        """
        self.H = np.random.random((self._num_bases, self._num_samples)) + 10**-4
        
    def _update_h(self):
        """ Overwrite for updating H.
        """
        pass

    def _update_w(self):
        """ Overwrite for updating W.
        """
        pass

    def _converged(self, i):
        """ 
        If the optimization of the approximation is below the machine precision,
        return True.

        Parameters
        ----------
            i   : index of the update step

        Returns
        -------
            converged : boolean 
        """
        derr = np.abs(self.ferr[i] - self.ferr[i-1])/self._num_samples
        if derr < self._EPS:
            return True
        else:
            return False

    def factorize(self, niter=100, show_progress=False, 
                  compute_w=True, compute_h=True, compute_err=True):
        """ Factorize s.t. WH = data
        
        Parameters
        ----------
        niter : int
                number of iterations.
        show_progress : bool
                print some extra information to stdout.
        compute_h : bool
                iteratively update values for H.
        compute_w : bool
                iteratively update values for W.
        compute_err : bool
                compute Frobenius norm |data-WH| after each update and store
                it to .ferr[k].
        
        Updated Values
        --------------
        .W : updated values for W.
        .H : updated values for H.
        .ferr : Frobenius norm |data-WH| for each iteration.
        """
        
        if show_progress:
            self._logger.setLevel(logging.INFO)
        else:
            self._logger.setLevel(logging.ERROR)        
        
        # create W and H if they don't already exist
        # -> any custom initialization to W,H should be done before
        if not hasattr(self,'W') and compute_w:
            self._init_w()
               
        if not hasattr(self,'H') and compute_h:
            self._init_h()                   
        
        # Computation of the error can take quite long for large matrices,
        # thus we make it optional.
        if compute_err:
            self.ferr = np.zeros(niter)
             
        for i in range(niter):
            if compute_w:
                self._update_w()

            if compute_h:
                self._update_h()                                        
         
            if compute_err:                 
                self.ferr[i] = self.frobenius_norm()                
                self._logger.info('FN: %s (%s/%s)'  %(self.ferr[i], i+1, niter))
            else:                
                self._logger.info('Iteration: (%s/%s)'  %(i+1, niter))
           

            # check if the err is not changing anymore
            if i > 1 and compute_err:
                if self._converged(i):
                    # adjust the error measure
                    self.ferr = self.ferr[:i]
                    break
                

class Kmeans(PyMFBase):
    """      
    Kmeans(data, num_bases=4)
    
    K-means clustering. Factorize a data matrix into two matrices s.t.
    F = | data - W*H | is minimal. H is restricted to unary vectors, W
    is simply the mean over the corresponding samples in "data".
    
    Parameters
    ----------
    data : array_like, shape (_data_dimension, _num_samples)
        the input data
    num_bases: int, optional
        Number of bases to compute (column rank of W and row rank of H).
        4 (default)     
    
    Attributes
    ----------
    W : "data_dimension x num_bases" matrix of basis vectors
    H : "num bases x num_samples" matrix of coefficients
    ferr : frobenius norm (after calling .factorize()) 
    
    Example
    -------
    Applying K-means to some rather stupid data set:
    
    >>> import numpy as np
    >>> data = np.array([[1.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    >>> kmeans_mdl = Kmeans(data, num_bases=2)
    >>> kmeans_mdl.factorize(niter=10)
    
    The basis vectors are now stored in kmeans_mdl.W, the coefficients in kmeans_mdl.H. 
    To compute coefficients for an existing set of basis vectors simply    copy W 
    to kmeans_mdl.W, and set compute_w to False:
    
    >>> data = np.array([[1.5], [1.2]])
    >>> W = np.array([[1.0, 0.0], [0.0, 1.0]])
    >>> kmeans_mdl = Kmeans(data, num_bases=2)
    >>> kmeans_mdl.W = W
    >>> kmeans_mdl.factorize(niter=1, compute_w=False)
    
    The result is a set of coefficients kmeans_mdl.H, s.t. data = W * kmeans_mdl.H.
    """        
    def _init_h(self):
        # W has to be present for H to be initialized  
        self.H = np.zeros((self._num_bases, self._num_samples))
        self._update_h()
         
    def _init_w(self):
        # set W to some random data samples
        sel = random.sample(range(self._num_samples), self._num_bases)
        
        # sort indices, otherwise h5py won't work
        self.W = self.data[:, np.sort(sel)]        
       
        
    def _update_h(self):                    
        # and assign samples to the best matching centers
        self.assigned = dist.vq(self.W, self.data)
        self.H = np.zeros(self.H.shape)
        self.H[self.assigned, range(self._num_samples)] = 1.0
                
                    
    def _update_w(self):            
        for i in range(self._num_bases):
            # cast to bool to use H as an index for data
            idx = np.array(self.H[i,:], dtype=np.bool)
            n = np.sum(idx)
            if n > 0:
                self.W[:,i] = np.sum(self.data[:, idx], axis=1)/n

class CNMF(PyMFBase):
    """      
    CNMF(data, num_bases=4)
    
    
    Convex NMF. Factorize a data matrix into two matrices s.t.
    F = | data - W*H | = | data - data*beta*H| is minimal. H and beta
    are restricted to convexity (beta >=0, sum(beta, axis=1) = [1 .. 1]).    
    
    Parameters
    ----------
    data : array_like, shape (_data_dimension, _num_samples)
        the input data
    num_bases: int, optional
        Number of bases to compute (column rank of W and row rank of H).
        4 (default)          
    
    Attributes
    ----------
    W : "data_dimension x num_bases" matrix of basis vectors
    H : "num bases x num_samples" matrix of coefficients
    ferr : frobenius norm (after calling .factorize()) 
    
    Example
    -------
    Applying CNMF to some rather stupid data set:
    
    >>> import numpy as np
    >>> from cnmf import CNMF
    >>> data = np.array([[1.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    >>> cnmf_mdl = CNMF(data, num_bases=2)
    >>> cnmf_mdl.factorize(niter=10)
    
    The basis vectors are now stored in cnmf_mdl.W, the coefficients in cnmf_mdl.H. 
    To compute coefficients for an existing set of basis vectors simply    copy W 
    to cnmf_mdl.W, and set compute_w to False:
    
    >>> data = np.array([[1.5, 1.3], [1.2, 0.3]])
    >>> W = [[1.0, 0.0], [0.0, 1.0]]
    >>> cnmf_mdl = CNMF(data, num_bases=2)
    >>> cnmf_mdl.W = W
    >>> cnmf_mdl.factorize(compute_w=False, niter=1)
    
    The result is a set of coefficients acnmf_mdl.H, s.t. data = W * cnmf_mdl.H.
    """

    # see .factorize() for the update of W and H
    # -> proper decoupling of W/H not possible ...
      
    def _init_h(self):
        if not hasattr(self, 'H'):
            # init basic matrices       
            self.H = np.zeros((self._num_bases, self._num_samples))
            
            # initialize using k-means
            km = Kmeans(self.data[:,:], num_bases=self._num_bases)        
            km.factorize(niter=10)
            assign = km.assigned
    
            num_i = np.zeros(self._num_bases)
            for i in range(self._num_bases):
                num_i[i] = len(np.where(assign == i)[0])
    
            self.H.T[range(len(assign)), assign] = 1.0                
            self.H += 0.2*np.ones((self._num_bases, self._num_samples))
        
        if not hasattr(self, 'G'):
            self.G = np.zeros((self._num_samples, self._num_bases))
            
            self.G[range(len(assign)), assign] = 1.0 
            self.G += 0.01                        
            self.G /= np.tile(np.reshape(num_i[assign],(-1,1)), self.G.shape[1])
            
        if not hasattr(self,'W'):
            self.W = np.dot(self.data[:,:], self.G)
    
    def _init_w(self):
        pass
    
    def factorize(self, niter=10, compute_w=True, compute_h=True, 
                  compute_err=True, show_progress=False):
        """ Factorize s.t. WH = data. 
            
            Parameters
            ----------
            niter : int
                    number of iterations.
            show_progress : bool
                    print some extra information to stdout.
            compute_h : bool
                    iteratively update values for H.
            compute_w : bool
                    iteratively update values for W.
            compute_err : bool
                    compute Frobenius norm |data-WH| after each update and store
                    it to .ferr[k].
            
            Updated Values
            --------------
            .W : updated values for W.
            .H : updated values for H.
            .ferr : Frobenius norm |data-WH| for each iteration.
        """   
        
        if not hasattr(self,'W'):
               self._init_w()
               
        if not hasattr(self,'H'):
                self._init_h()              
        
        def separate_positive(m):
            return (np.abs(m) + m)/2.0 
        
        def separate_negative(m):
            return (np.abs(m) - m)/2.0 
            
        if show_progress:
            self._logger.setLevel(logging.INFO)
        else:
            self._logger.setLevel(logging.ERROR)                 
            
        XtX = np.dot(self.data[:,:].T, self.data[:,:])
        XtX_pos = separate_positive(XtX)
        XtX_neg = separate_negative(XtX)
        
        self.ferr = np.zeros(niter)
        # iterate over W and H
        
        for i in range(niter):
            # update H
            XtX_neg_x_W = np.dot(XtX_neg, self.G)
            XtX_pos_x_W = np.dot(XtX_pos, self.G)
            
            if compute_h:
                H_x_WT = np.dot(self.H.T, self.G.T)                                
                ha = XtX_pos_x_W + np.dot(H_x_WT, XtX_neg_x_W)
                hb = XtX_neg_x_W + np.dot(H_x_WT, XtX_pos_x_W) + 10**-9
                self.H = (self.H.T*np.sqrt(ha/hb)).T
        
            # update W            
            if compute_w:
                HT_x_H = np.dot(self.H, self.H.T)
                wa = np.dot(XtX_pos, self.H.T) + np.dot(XtX_neg_x_W, HT_x_H)
                wb = np.dot(XtX_neg, self.H.T) + np.dot(XtX_pos_x_W, HT_x_H) + 10**-9
            
                self.G *= np.sqrt(wa/wb)
                self.W = np.dot(self.data[:,:], self.G)
                                
            if compute_err:                 
                self.ferr[i] = self.frobenius_norm()
                self._logger.info('FN: %s (%s/%s)'  %(self.ferr[i], i+1, niter))
            else:                
                self._logger.info('Iteration: (%s/%s)'  %(i+1, niter))

            if i > 1 and compute_err:
                if self._converged(i):                    
                    self.ferr = self.ferr[:i]                    
                    break


#-------------------------------------------------------------------------------
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
        model = CNMF(np.transpose(data), num_bases=k, niter=niter)
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

