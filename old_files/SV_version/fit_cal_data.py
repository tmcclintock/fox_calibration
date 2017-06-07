"""
Here I fit the calibration data.

Note that the calibration data is stored locally, since the covariance
matrix hasn't been released yet by DES, so it cannot be put on github.
"""
import numpy as np
import os, sys, emcee
import scipy.optimize as op
import matplotlib.pyplot as plt
sys.path.insert(0, "../../../Delta-Sigma/src/wrapper/")
import py_Delta_Sigma as pyDS
from colossus.halo import concentration as conc
from colossus.cosmology import cosmology as col_cosmology

Y1 = False

#Paths to the data and covariance
if Y1:
    print "not implemented yet"
else:
    datapath = "/home/tmcclintock/Desktop/des_wl_work/SV_work/calibration_work/calibration_data/cal_profile_z%.2f_m%d.txt"
    covpath  = "/home/tmcclintock/Desktop/des_wl_work/SV_work/calibration_work/calibration_data/cal_cov_z%.2f_m%d.txt"

true_mass = np.loadtxt("txt_files/SV_masses.txt")#True masses
true_lM = np.log10(true_mass)

def get_richind():
    """
    Get the index of the covariance matrix that we will use for each mass bin.
    They don't line up exactly, like how the Y1 bins are designed.
    """
    lam_edges = [5, 10, 14, 20, 30, 45, 60, np.inf]

    lM = np.log10(true_mass)
    rich = 30.*10**((lM-14.37)/1.12) #approximate richnesses
    inds = np.zeros_like(rich)
    for i in range(len(inds)):
        for j in range(len(inds[i])):
            for k in range(len(lam_edges)-1):
                if rich[i,j] < lam_edges[k+1]:
                    inds[i,j] = k
                    break
    return inds

#Snapshot characteristics
zs = [1.0, 0.5, 0.25, 0.0]
if Y1:
    inds = range(0, 7)
else:
    edges = [13.1, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.5, 15.0]#, 16.0]
    inds = get_richind()
    print inds
Nz = len(zs)
Nl = len(inds[0])

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
                "miscentering":0,"averaging":1,
                "single_miscentering":0}

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
    return -0.5*np.dot(X, np.dot(icov, X))

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
    for i in range(len(zs)):
        z = zs[i]
        #Get power spectra
        k = np.loadtxt("../txt_files/P_files/k.txt")
        Plin = np.loadtxt("../txt_files/P_files/Plin_z%.2f.txt"%z)
        Pmm  = np.loadtxt("../txt_files/P_files/Pnl_z%.2f.txt"%z)
        extras = (k, Plin, Pmm, cosmo, input_params)
        #Give the bin edges in comoving units; Mpc/h
        input_params["R_bin_min"] = 0.0323*(h*(1+z))
        input_params["R_bin_max"] = 30.0*(h*(1+z))
        for j in range(len(inds[0])):
            covindex = inds[i,j]
            #Data is in units of Msun/pc^2 physical
            R, DS, err, flag = np.loadtxt(datapath%(z, j), unpack=True)
            flags = np.where(flag==1)[0]
            cov = np.loadtxt(covpath%(z, covindex))
            cov = cov[flags]
            cov = cov[:, flags]
            icov = np.linalg.inv(cov)            
            nll = lambda *args: -lnprob(*args)
            result = op.minimize(nll, x0=true_lM[i,j],
                                 args=(R, DS, icov, flags, z, extras))
            print result
            bfmasses[i, j] = result['x'][0]
    return bfmasses

def do_mcmc(bfmasses):
    """
    Find the mcmc masses, where we will end up starting the chain.
    """
    nwalkers = 4
    ndim = 1
    nsteps = 1000
    for i in range(len(zs)):
        z = zs[i]
        #Get power spectra
        k = np.loadtxt("../txt_files/P_files/k.txt")
        Plin = np.loadtxt("../txt_files/P_files/Plin_z%.2f.txt"%z)
        Pmm  = np.loadtxt("../txt_files/P_files/Pnl_z%.2f.txt"%z)
        #Give the bin edges in comoving units; Mpc/h
        input_params["R_bin_min"] = 0.0323*(h*(1+z))
        input_params["R_bin_max"] = 30.0*(h*(1+z))
        extras = (k, Plin, Pmm, cosmo, input_params)
        for j in range(len(inds[0])):
            covindex = inds[i,j]
            pos = [bfmasses[i, j] + 1e-4*np.random.randn(ndim) for nw in range(nwalkers)]
            R, DS, err, flag = np.loadtxt(datapath%(z, j), unpack=True)
            flags = np.where(flag==1)[0]
            cov = np.loadtxt(covpath%(z, covindex))
            cov = cov[flags]
            cov = cov[:, flags]
            icov = np.linalg.inv(cov)
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(R, DS, icov, flags, z, extras), threads=8)
            print "Running chain for z=%f m%d starting at lM=%.3f"%(z, j, bfmasses[i,j])
            sampler.run_mcmc(pos, nsteps)
            print "\tchain complete for z=%f m%d starting at lM=%.3f"%(z, j, bfmasses[i,j])
            chain = sampler.flatchain
            likes = sampler.flatlnprobability
            np.savetxt("txt_files/chains/chain_z%.2f_m%d.txt"%(z, j), chain)
            np.savetxt("txt_files/chains/likes_z%.2f_m%d.txt"%(z, j), likes)
    return 0

def reduce_chains():
    nsteps = 100
    nburn  = nsteps/2
    nwalkers = 4
    masses = np.ones((Nz,Nl))
    err  = np.ones_like(masses)
    for i in range(len(zs)):
        z = zs[i]
        for j in inds:
            chain = np.genfromtxt("txt_files/chains/chain_z%.2f_m%d.txt"%(z, j))
            masses[i, j] = np.mean(chain)
            err[i, j]  = np.std(chain)
            print "lM=%.3f +- %.3f at z%.2f m%d"%(masses[i,j], err[i,j], z, j)
    np.savetxt("txt_files/SV_mcmc_masses.txt", masses)
    np.savetxt("txt_files/SV_mcmc_errs.txt", err)
    return

if __name__ == "__main__":
    #bfmasses = best_fits()
    #np.savetxt("txt_files/SV_BF_masses.txt", bfmasses)

    bfmasses = np.loadtxt("txt_files/SV_BF_masses.txt")
    do_mcmc(bfmasses)

    #reduce_chains()
