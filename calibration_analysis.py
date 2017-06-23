"""
This is a file that combines all the small scripts I used earlier.
This is being written so that I can use general paths to the catalogs of 
interest isntead of hard-coding in paths to certain catalogs.

Note: the 'richness calibration' is in old_files/ and the 'mass calibration'
is in old_files/SV_version.
"""
import os, sys
import numpy as np
sys.path.insert(0, "./src/")
from CF_functions import *
from DS_functions import *

inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
halobase="/calvin1/tmcclintock/fox_data/richness_halos/rich_snapdir_ps%d_%03d/"
mpath = halobase+"/mass_halos_ps%d_m%d_%03d.txt" #Path to mass split halos
lpath = halobase+"/richness_halos_ps%d_l%d_%03d.txt" #Path to lam split halos
lam_edges = [5, 10, 14, 20, 30, 45, 60, np.inf]
lM_edges = [13.1, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.5, 15.0]#, 16.0]
mmeans = np.genfromtxt("L_ps25_masses.txt")

dmpath = "/calvin1/tmcclintock/down_sampled_snapshots/snapdir_%03d/snapshot_%03d_z%s_down10000"

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    """
    for ps in [15, 25, 35]:
        for i,ind in zip(range(len(inds)), inds):
            z = zs[i]
            for j in linds:
                halos = lpath%(ps, ind, ps, j, ind)
                calc_hhcf(halos, "output_files/hhcf/hhcf_ps%d_z%d_l%d.txt"%(ps,i,j))
                print "HHCF ps%d z%d, l%d complete"%(ps,i,j)
                continue
            continue
        continue
    print "HHCFs created"
    """
    #RR, DdRd, Dd, Rh = calc_DdRd_and_RR(dmpath%(6,6,"1.0"))
    print "DdRd and RR used in HMCF done"
    #R,xihm = calc_hmcf(lpath%(pscatter,6,pscatter,3,6), RR, DdRd, Dd, Rh, "test.txt")
    R, xihm = np.loadtxt("test.txt")
    inds = np.invert(np.isnan(xihm))
    R = R[inds]
    xihm = xihm[inds]
    print "HMCF created"

    #Now with the known mean masses, call Build_Delta_Sigma
    Mass = mmeans[0,3]
    redshift = 1.0
    R, DS = calc_DS(R, xihm, Mass, redshift, "dstext.txt")
    print "DeltaSigma built"

    #Now pair it up with a covariance matrix and output the correct files
    #Read in the associated covariance matrix
    create_data_vector(R, DS, cov, z, save, Csave)
    print "Data vector created"
