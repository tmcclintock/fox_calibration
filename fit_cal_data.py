"""
Here I fit the calibration data.

Note that the calibration data is stored locally, since the covariance
matrix hasn't been released yet by DES, so it cannot be put on github.
"""
import numpy as np
import os, sys, emcee
import scipy.optimize as op
import matplotlib.pyplot as plt
sys.path.insert(0, "../Delta-Sigma/src/wrapper/")
import py_Delta_Sigma as pyDS

#Paths to the data and covariance
datapath = "/home/tmcclintock/Desktop/des_wl_work/Y1_work/blinded_data/calibration_data/cal_profile_z%.2f_l%d.txt"
covpath = "/home/tmcclintock/Desktop/des_wl_work/Y1_work/blinded_data/calibration_data/cal_cov_z%.2f_l%d.txt"

#Snapshot characteristics
inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
linds = range(0,7)
Nz = len(zs)
Nl = len(linds)

#Get the power spectra from somewhere TODO

#Write the likelihoods
def lnprior(lM):
    if lM < 10: return -np.inf
    return 0.0

def lnlike(lM, R, ds, icov):
    model = np.ones_like(R) #replace
    X = ds - model
    return np.dot(X, np.dot(icov, X))

def lnprob(lM, R, ds, icov):
    lpr = lnprior(lM)
    if not np.isfinite(lpr): return -np.inf
    return lpr + lnlike(lM, R, ds, icov)

#Define the functions that we use for flow control
def best_fits():
    """
    Find the best fit masses, where we will end up starting the chain.
    """
    bfmasses = np.ones((Nz,Nl))
    for i in range(len(inds)):
        index = inds[i]
        z = zs[i]
        for j in linds:
            R, DS, err, flag = np.loadtxt(datapath%(z, j), unpack=True)
            print R.shape
            sys.exit()
    return bfmasses

if __name__ == "__main__":
    best_fits()
