"""
Plot the best fit and mcmc results.
"""
import numpy as np
import os, sys, emcee
import scipy.optimize as op
import matplotlib.pyplot as plt
plt.rc("text", usetex=True, fontsize=24)
sys.path.insert(0, "../Delta-Sigma/src/wrapper/")
import py_Delta_Sigma as pyDS
from colossus.halo import concentration as conc
from colossus.cosmology import cosmology as col_cosmology

#Paths to the data and covariance
datapath = "/home/tmcclintock/Desktop/des_wl_work/Y1_work/blinded_data/calibration_data/cal_profile_z%.2f_l%d.txt"

#Snapshot characteristics
inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
linds = range(0,7)
Nl = len(linds)
c = np.linspace(1.0, 0.3, Nl)
cmaps = ['Reds','Oranges','Greens', 'Blues']

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
params = {"NR":300,"Rmin":0.01,
          "Rmax":200.0,"Nbins":15,"delta":200,
          "Rmis":0.2, "fmis":0.0,
          "miscentering":0,"averaging":1}

#masses = np.loadtxt("txt_files/BF_masses.txt")
masses = np.loadtxt("txt_files/mcmc_masses.txt")
#masses = np.log10(np.loadtxt("txt_files/true_masses.txt"))

for i in range(len(inds)):
    cmap = plt.get_cmap(cmaps[i])
    index = inds[i]
    z = zs[i]
    #Get power spectra
    Plin = np.loadtxt("txt_files/P_files/Plin_z%.2f.txt"%z)
    Pnl  = np.loadtxt("txt_files/P_files/Pnl_z%.2f.txt"%z)
    #Give the bin edges in comoving units; Mpc/h
    params["R_bin_min"] = 0.0323*h*(1+z)
    params["R_bin_max"] = 30.0*h*(1+z)
    for j in linds:
        R, DS, err, flag = np.loadtxt(datapath%(z, j), unpack=True)
        lM = masses[i,j]
        print 10**lM
        params['Mass'] = 10**lM
        params["concentration"] = conc.concentration(10**lM, 
                                                     '200m', z, 
                                                     model='diemer15')
        result = pyDS.calc_Delta_Sigma(k, Plin, k, Pnl, cosmo, params)
        model = result['ave_delta_sigma']*h*(1.+z)**2 #Msun/pc^2 physical
        plt.loglog(R, model, c=cmap(c[j]), ls='-')
        plt.errorbar(R, DS, err, c=cmap(c[j]), marker='o', ls='')
    plt.ylim(.1, 1e3)
    plt.xlabel(r"$R\ [{\rm Mpc}]$")
    plt.ylabel(r"$\Delta\Sigma\ [{\rm M_\odot/pc^2}]$")
    plt.title("z=%.2f"%z)
    plt.subplots_adjust(bottom=0.17, left=0.2)
    plt.gcf().savefig("BFs_z%.2f.png"%z)
    plt.show()
