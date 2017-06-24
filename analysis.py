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

inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
DSdatabase      = "/calvin1/tmcclintock/DES_Y1_data/calibration_data/ds_ps%d_z%d_l%d.txt"
covdatabase     = "/calvin1/tmcclintock/DES_Y1_data/calibration_data/cov_ps%d_z%d_l%d.txt"

#Fox cosmology
h = 0.670435
cosmo = {"h":h,"om":0.31834,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]
input_params = {"NR":300,"Rmin":0.01,
                "Rmax":200.0,"Nbins":15,"delta":200,
                "Rmis":0.0, "fmis":0.0,
                "miscentering":0,"averaging":1}
k = np.loadtxt("txt_files/P_files/k.txt")

if __name__ == "__main__":
    print "not implemented"

    for ps in [15]:#, 25, 35]:
        mean_masses = np.genfromtxt("L_ps%d_masses.txt"%ps)
        bf_masses = np.zeros_like(mean_masses)
        for i,ind in zip(range(len(inds)), inds):
            z = zs[i]
            Plin = np.loadtxt("txt_files/P_files/Plin_z%.2f.txt"%z)
            Pmm  = np.loadtxt("txt_files/P_files/Pnl_z%.2f.txt"%z)
            extras = [k, Plin, Pmm, cosmo, input_params]
            for j in linds:
                DSpath  = DSdatabase%(ps, i, j)
                covpath = covdatabase%(ps, i, j)
                """
                R, DS = np.loadtxt(DSpath)
                cov = np.loadtxt(covpath)
                inds = R<0.2 #Mpc
                cov = cov[inds]
                cov = cov[:,inds]
                icov = np.linalg.inv(cov)
                DS = DS[inds]
                R = R[inds]
                """
                """
                nll = lambda *args: -lnprob(*args)
                result = op.minimize(nll, x0=true_lM[i,j],args=(R, DS, icov, inds, z, extras))
                bf_masses[i,j] = result['x']
                """
                print "Best fit done for ps%d z%d, l%d"%(ps, ind, j)
                continue #end j
            #np.savetxt("output_files/mass_fits/bf_masses_ps%d.txt", bf_masses)
            continue #end i,ind
        continue #end ps
