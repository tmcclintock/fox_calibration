"""
Here I actually do the calibration analysis, calling the optimizer/MCMC for each redshift, richness, and scatter value.
"""
import os, sys
sys.path.insert(0,"./src/")
from likelihoods import *
import numpy as np
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op

#MCMC parameters
nwalkers = 4
ndim = 1
nsteps = 4000
nburn = nsteps/4

inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
DSdatabase      = "/calvin1/tmcclintock/DES_Y1_data/calibration_data/dss/ds_ps%d_z%d_l%d.txt"
covdatabase     = "/calvin1/tmcclintock/DES_Y1_data/calibration_data/covs/cov_ps%d_z%d_l%d.txt"
svcovdatabase   = "/calvin1/tmcclintock/DES_Y1_data/calibration_data/svcovs/svcov_ps%d_z%d_l%d.txt"
use_y1 = True

#Fox cosmology
h = 0.670435
cosmo = {"h":h,"om":0.31834,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]
input_params = {"NR":300,"Rmin":0.01,
                "Rmax":200.0,"Nbins":15,"delta":200,
                "Rmis":0.0, "fmis":0.0,
                "miscentering":0,"averaging":1,"single_miscentering":0}
use_old_P = False #We can choose to use non-Takahashi P_mm by making this True
klin = np.loadtxt("txt_files/P_files/k.txt")
if use_old_P: knl  = np.loadtxt("txt_files/P_files/knl.txt")
else: knl = np.loadtxt("txt_files/P_files/k.txt")

def do_best_fit(ps, zi, lj):
    true_lM = np.log10(np.genfromtxt("L_ps%d_masses.txt"%ps))[zi, lj]
    z = zs[zi]
    Plin = np.loadtxt("txt_files/P_files/Plin_z%.2f.txt"%z)
    if use_old_P: Pmm  = np.loadtxt("txt_files/P_files/Pnl_old_z%.2f.txt"%z)
    else: Pmm  = np.loadtxt("txt_files/P_files/Pnl_z%.2f.txt"%z)
    input_params["R_bin_min"] = 0.0323*(h*(1+z)) #Mpc/h comoving
    input_params["R_bin_max"] = 30.0*(h*(1+z)) #Mpc/h comoving
    extras = [klin, knl, Plin, Pmm, cosmo, input_params]
    DSpath  = DSdatabase%(ps, zi, lj)
    if use_y1: covpath = covdatabase%(ps, zi, lj)
    else: covpath = svcovdatabase%(ps, zi, lj)
    R, DS = np.loadtxt(DSpath).T
    cov = np.loadtxt(covpath)
    cut = R>0.2 #Mpc
    cov = cov[cut]
    cov = cov[:,cut]
    icov = np.linalg.inv(cov)
    DS = DS[cut]
    R = R[cut]
    nll = lambda *args: -lnprob(*args)
    result = op.minimize(nll, x0=true_lM,args=(R, DS, icov, cut, z, extras))
    bf_lM = result['x']
    bf_cal = 10**true_lM/10**bf_lM
    if use_y1:  print "Best fit done for ps%d z%d, l%d, C_y1"%(ps, inds[zi], lj)
    else:  print "Best fit done for ps%d z%d, l%d, C_sv"%(ps, inds[zi], lj)
    print "Bf is ",bf_lM, bf_cal
    return 10**bf_lM
"""
    for ps in [15, 25, 35]:
        true_lM = np.log10(np.genfromtxt("L_ps%d_masses.txt"%ps))
        bf_masses = np.zeros_like(true_lM)
        cal = np.zeros_like(bf_masses)
        for i,ind in zip(range(len(inds)), inds):
            z = zs[i]
            Plin = np.loadtxt("txt_files/P_files/Plin_z%.2f.txt"%z)
            if use_old_P: Pmm  = np.loadtxt("txt_files/P_files/Pnl_old_z%.2f.txt"%z)
            else: Pmm  = np.loadtxt("txt_files/P_files/Pnl_z%.2f.txt"%z)
            input_params["R_bin_min"] = 0.0323*(h*(1+z)) #Mpc/h comoving
            input_params["R_bin_max"] = 30.0*(h*(1+z)) #Mpc/h comoving
            extras = [klin, knl, Plin, Pmm, cosmo, input_params]
            for j in linds:
                DSpath  = DSdatabase%(ps, i, j)
                covpath = covdatabase%(ps, i, j)
                R, DS = np.loadtxt(DSpath).T
                cov = np.loadtxt(covpath)
                cut = R>0.2 #Mpc
                cov = cov[cut]
                cov = cov[:,cut]
                icov = np.linalg.inv(cov)
                DS = DS[cut]
                R = R[cut]
                nll = lambda *args: -lnprob(*args)
                result = op.minimize(nll, x0=true_lM[i,j],args=(R, DS, icov, cut, z, extras))

                bf_masses[i,j] = result['x']
                cal[i,j] = 10**true_lM[i,j]/10**bf_masses[i,j]
                print "Best fit done for ps%d z%d, l%d"%(ps, ind, j)
                print "Bf is ",result['x'], cal[i,j]
                continue #end j
            continue #end i,ind
        np.savetxt("output_files/mass_fits/bf_masses_ps%d.txt"%ps, bf_masses)
        np.savetxt("output_files/mass_fits/bf_cal_ps%d.txt"%ps, cal)
        continue #end ps
"""

def do_mcmc():
    for ps in [35]:#[15, 25, 35]:
        true_lM = np.log10(np.genfromtxt("L_ps%d_masses.txt"%ps))
        bf_masses = np.loadtxt("output_files/mass_fits/bf_masses_ps%d.txt"%ps)
        mcmc_masses = np.zeros_like(true_lM)
        mcmc_stds   = np.zeros_like(mcmc_masses)
        cal     = np.zeros_like(mcmc_masses)
        calerr = np.zeros_like(mcmc_masses)
        for i,ind in zip(range(len(inds)), inds):
            z = zs[i]
            Plin = np.loadtxt("txt_files/P_files/Plin_z%.2f.txt"%z)
            if use_old_P: Pmm  = np.loadtxt("txt_files/P_files/Pnl_old_z%.2f.txt"%z)
            else: Pmm  = np.loadtxt("txt_files/P_files/Pnl_z%.2f.txt"%z)
            input_params["R_bin_min"] = 0.0323*(h*(1+z))
            input_params["R_bin_max"] = 30.0*(h*(1+z))
            extras = [klin, knl, Plin, Pmm, cosmo, input_params]
            for j in linds:
                DSpath  = DSdatabase%(ps, i, j)
                covpath = covdatabase%(ps, i, j)
                R, DS = np.loadtxt(DSpath).T
                cov = np.loadtxt(covpath)
                cut = R>0.2 #Mpc
                cov = cov[cut]
                cov = cov[:,cut]
                icov = np.linalg.inv(cov)
                DS = DS[cut]
                R = R[cut]
                sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(R, DS, icov, cut, z, extras), threads=4)
                pos = [bf_masses[i, j] + 1e-4*np.random.randn(ndim) for nw in range(nwalkers)]
                print "Running chain for z=%f l%d starting at lM=%.3f"%(z, j, bf_masses[i,j])
                sampler.run_mcmc(pos, nsteps)
                print "\tchain complete for z=%f l%d starting at lM=%.3f"%(z, j, bf_masses[i,j])
                fullchain = sampler.flatchain
                likes = sampler.flatlnprobability
                np.savetxt("output_files/chains/chain_ps%d_z%.2f_l%d.txt"%(ps, z, j), fullchain)
                np.savetxt("output_files/chains/likes_ps%d_z%.2f_l%d.txt"%(ps, z, j), likes)
                chain = fullchain[nburn*nwalkers:]
                mcmc_masses[i,j] = np.mean(chain)
                mcmc_stds[i,j] = np.std(chain)
                cal[i,j] = 10**true_lM[i,j]/10**mcmc_masses[i,j]
                calerr[i,j] = np.log(10) * mcmc_stds[i,j] * cal[i,j]
                print "MCMC done for ps%d z%d, l%d"%(ps, i, j)
                print "mcmc result is ",mcmc_masses[i,j], cal[i,j]
                continue #end j
            continue #end i,ind
        np.savetxt("output_files/mass_fits/mcmc_masses_ps%d.txt"%ps, mcmc_masses)
        np.savetxt("output_files/mass_fits/mcmc_cal_ps%d.txt"%ps, cal)
        np.savetxt("output_files/mass_fits/mcmc_calerr_ps%d.txt"%ps, calerr)
        continue #end ps

if __name__ == "__main__":
    save = False
    for ps in [0]:#[15, 25, 35]:
        true_M = np.genfromtxt("L_ps%d_masses.txt"%ps)
        bf_M   = np.ones_like(true_M)
        bf_cal = np.ones_like(true_M)
        for i in xrange(3,4):#len(inds)):
            for j in xrange(3,7):#linds:
                use_y1 = True
                bf_M[i, j] = do_best_fit(ps, i, j)
                use_y1 = False
                bf_M[i, j] = do_best_fit(ps, i, j)
                bf_cal[i, j] = true_M[i, j]/bf_M[i, j]
        if save: 
            np.savetxt("output_files/mass_fits/bf_masses_ps%d.txt"%ps, bf_M)
            np.savetxt("output_files/mass_fits/bf_cal_ps%d.txt"%ps, bf_cal)
    #do_mcmc()
    #see_results() #NEED TO IMPLEMENT
