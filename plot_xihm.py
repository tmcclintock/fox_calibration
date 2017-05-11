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

#The radial points of all data
nbins = 50
bins = np.logspace(np.log10(0.01), np.log10(150.0), nbins+1)
Rd = (bins[:-1]+bins[1:])/2.

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

masses = np.loadtxt("txt_files/true_masses.txt")

for i in range(len(inds)):
    if i < 0: continue
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
        lM = np.log10(masses[i,j])
        print 10**lM
        params['Mass'] = 10**lM
        params["concentration"] = conc.concentration(10**lM, 
                                                     '200m', z, 
                                                     model='diemer15')
        k2 = np.copy(k)
        result = pyDS.calc_Delta_Sigma(k, Plin, k2, Pnl, cosmo, params)
        R = result['R']
        xi = result['xi_hm']
        xi_data = np.loadtxt("txt_files/richness_txt_files/hmcf_z%.2f_l%d.txt"%(z, j))
        #print xi.shape, xi_data.shape, R.shape, Rd.shape
        plt.loglog(R, xi, label=r"Model", c='r')
        #plt.loglog(R, result['xi_1halo'], c='b', ls='--')
        #plt.loglog(R, result['xi_2halo'], c='b', ls='--')
        plt.loglog(Rd, xi_data, label=r"Sims", c='g')
        plt.xlabel(r"$R\ [{\rm Mpc/h}]$")
        plt.ylabel(r"$\xi_{\rm hm}$")
        plt.subplots_adjust(bottom=0.15, left=0.15)
        #plt.legend()
        plt.title("z=%.2f"%z)
        plt.gcf().savefig("xi_comparison_z%.2f.png"%z)
    plt.show()
