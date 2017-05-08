"""
Use the Build-Delta-Sigma repository to turn my HMCF curves into delta sigma profiles.
"""
import os, sys
sys.path.insert(0, "../Build-Delta-Sigma/src/wrapper/")
import py_Build_Delta_Sigma
import matplotlib.pyplot as plt
import numpy as np
from colossus.halo import concentration as conc
from colossus.cosmology import cosmology as col_cosmology

plt.rc("text",usetex=True, fontsize=24)

inpath  = "txt_files/richness_txt_files/hmcf_z%0.2f_l%d.txt"

#The radial points of all data
nbins = 50
bins = np.logspace(np.log10(0.01), np.log10(150.0), nbins+1)
R = (bins[:-1]+bins[1:])/2.

#Snapshot characteristics
inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
c = np.linspace(0.5, 0.2, len(linds))
cmaps = ['Reds','Oranges','Greens', 'Blues']

#This is the fox sim cosmology
cosmo = {"h":0.670435,"om":0.31834,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]
#Here is the same cosmology but for cosmocalc
colcos = {"H0":cosmo['h']*100.,"Om0":cosmo['om'], 
          'Ob0': 0.049017, 'sigma8': 0.83495, 'ns': 0.96191, 'flat':True}
col_cosmology.addCosmology('fiducial_cosmology', colcos)
col_cosmology.setCosmology('fiducial_cosmology')

#Get the mean masses
masses = np.loadtxt("txt_files/mean_masses.txt")

#Loop over all bins
DeltaSigmas = []
for i in range(len(inds)):
    index = inds[i]
    z = zs[i]
    zstring = zstrings[i]
    cmap = plt.get_cmap(cmaps[i])
    for j in linds:
        M = masses[i,j]
        xi_hm = np.loadtxt(inpath%(z, j))
        input_params = {"Mass": masses[i, j], "delta":200, 
                        "timing":1, "miscentering":0}
        input_params["concentration"] = conc.concentration(M, '200m', z, 
                                                           model='diemer15')

        results = py_Build_Delta_Sigma.build_Delta_Sigma(R, xi_hm, 
                                                          cosmo, 
                                                          input_params)
        ds = results['delta_sigma']
        print results.keys()
        if j == 0:
            plt.loglog(R, ds, c=cmap(c[j]), label=r"z=%.2f"%(z))
        else:
            plt.loglog(R, ds, c=cmap(c[j]))#, label=r"z=%.2f l%d"%(z,j))
        np.savetxt("txt_files/richness_txt_files/deltasigma_z%.2f_l%d.txt"%(z, j), ds)
    plt.legend(loc=0, fontsize=12)
plt.xlabel(r"$R\ [{\rm Mpc}/h]$")
plt.ylabel(r"$\Delta\Sigma\ [{\rm M_\odot}h/{\rm pc^2}]$")
plt.subplots_adjust(bottom=0.15, left=0.15)
#plt.xlim(0.1, 1000)
plt.show()
