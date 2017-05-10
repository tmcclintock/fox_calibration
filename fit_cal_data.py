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
from colossus.halo import concentration as conc
from colossus.cosmology import cosmology as col_cosmology

#Paths to the data and covariance
datapath = "/home/tmcclintock/Desktop/des_wl_work/Y1_work/blinded_data/calibration_data/cal_profile_z%.2f_l%d.txt"
covpath = "/home/tmcclintock/Desktop/des_wl_work/Y1_work/blinded_data/calibration_data/cal_cov_z%.2f_l%d.txt"

#Snapshot characteristics
inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
linds = range(0,7)
Nz = len(zs)
Nl = len(linds)

#Get the power spectra wavenumbers
k = np.loadtxt("txt_files/P_files/k.txt")

#This is the fox sim cosmology
h = 0.670435
cosmo = {"h":h,"om":0.31834,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]
colcos = {"H0":cosmo['h']*100.,"Om0":cosmo['om'], 
          'Ob0': 0.049017, 'sigma8': 0.83495, 'ns': 0.96191, 'flat':True}
col_cosmology.addCosmology('fiducial_cosmology', colcos)
col_cosmology.setCosmology('fiducial_cosmology')
input_params = {"NR":300,"Rmin":0.01,
                "Rmax":200.0,"Nbins":15,"delta":200,
                "Rmis":0.0, "fmis":0.0,
                "miscentering":0,"averaging":1}

#Write the likelihoods
def lnprior(params):
    lM = params
    #Mass to low or too high
    if lM < 10 or lM > 18: return -np.inf
    return 0.0

def lnlike(params, R, ds, icov, flags, z, extras):
    lM = params
    k, Plin, Pnl, cosmo, inparams = extras
    inparams['Mass'] = 10**lM #Mpc/h
    inparams["concentration"] = conc.concentration(10**lM, '200m', z, model='diemer15')
    result = pyDS.calc_Delta_Sigma(k, Plin, k, Pnl, cosmo, inparams)
    model = result['ave_delta_sigma']*h*(1.+z)**2 #Msun/pc^2 physical
    X = ds - model
    X = X[flags]
    return np.dot(X, np.dot(icov, X))

def lnprob(params, R, ds, icov, flags, z, extras):
    lM = params
    lpr = lnprior(params)
    if not np.isfinite(lpr): return -np.inf
    return lpr + lnlike(params, R, ds, icov, flags, z, extras)

#Define the functions that we use for flow control
def best_fits():
    """
    Find the best fit masses, where we will end up starting the chain.
    """
    bfmasses = np.ones((Nz,Nl))
    for i in range(len(inds)):
        index = inds[i]
        z = zs[i]
        #Get power spectra
        Plin = np.loadtxt("txt_files/P_files/Plin_z%.2f.txt"%z)
        Pmm  = np.loadtxt("txt_files/P_files/Pnl_z%.2f.txt"%z)
        extras = (k, Plin, Pmm, cosmo, input_params)
        #Give the bin edges in comoving units; Mpc/h
        input_params["R_bin_min"] = 0.0323*(h*(1+z))
        input_params["R_bin_max"] = 30.0*(h*(1+z))
        for j in linds:
            R, DS, err, flag = np.loadtxt(datapath%(z, j), unpack=True)
            flags = np.where(flag==1)[0]
            cov = np.loadtxt(covpath%(z, j))
            cov = cov[flags]
            cov = cov[:, flags]
            icov = np.linalg.inv(cov)            
            nll = lambda *args:lnprob(*args)
            result = op.minimize(nll, x0=13.2,
                                 args=(R, DS, icov, flags, z, extras))
            print result
            bfmasses[i, j] = result['x'][0]
    return bfmasses

def do_mcmc():
    return 0

if __name__ == "__main__":
    bfmasses = best_fits()
    np.savetxt("txt_files/BF_masses_withc.txt", bfmasses)

    bfmasses = np.loadtxt("txt_files/BF_masses.txt")
    print bfmasses
    do_mcmc()
