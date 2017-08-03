"""
This is a file that combines all the small scripts I used earlier.
This is being written so that I can use general paths to the catalogs of 
interest isntead of hard-coding in paths to certain catalogs.

Note: the 'richness calibration' is in old_files/ and the 'mass calibration'
is in old_files/SV_version.
"""
import os, sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "./src/")
from CF_functions import *
from DS_functions import *

inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
halobase= "/calvin1/tmcclintock/fox_data/richness_halos/rich_snapdir_ps%d_%03d/"
covbase = "/calvin1/tmcclintock/DES_Y1_data/tamas_v1/full-mcal-urem_y1subtr_l%d_z%d_dst_cov.dat"
mpath = halobase+"/mass_halos_ps%d_m%d_%03d.txt" #Path to mass split halos
lpath = halobase+"/richness_halos_ps%d_l%d_%03d.txt" #Path to lam split halos
lam_edges = [5, 10, 14, 20, 30, 45, 60, np.inf]
lM_edges = [13.1, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.5, 15.0]#, 16.0]

dmpath_base = "/calvin1/tmcclintock/down_sampled_snapshots/snapdir_%03d/snapshot_%03d_z%s_down100"

hhcf_savepath = "output_files/hhcf/hhcf_ps%d_z%d_l%d.txt"
hmcf_savepath = "output_files/hmcf/hmcf_ps%d_z%d_l%d.txt"
ds_savepath   = "output_files/ds/ds_ps%d_z%d_l%d.txt"
DSdatabase      = "/calvin1/tmcclintock/DES_Y1_data/calibration_data/ds_ps%d_z%d_l%d.txt"
covdatabase     = "/calvin1/tmcclintock/DES_Y1_data/calibration_data/cov_ps%d_z%d_l%d.txt"
covinds = [2, 2, 2, 1] #Covariance matrices to use for each snap

if __name__ == "__main__":
    """
    for ps in [15, 25, 35]:
        for i,ind in zip(range(len(inds)), inds):
            z = zs[i]
            for j in linds:
                halos = lpath%(ps, ind, ps, j, ind)
                calc_hhcf(halos, hhcf_savepath%(ps,i,j))
                print "HHCF ps%d z%d, l%d complete"%(ps,ind,j)
                continue #end j
            continue #end i,ind
        continue #end ps
    print "HHCFs created"
    """
    
    for i,ind in zip(range(len(inds)), inds):
        z = zs[i]
        zstring = zstrings[i]
        dmpath = dmpath_base%(ind, ind, zstring)
        RR, DdRd, Dd, Rh = calc_DdRd_and_RR(dmpath)
        print "DdRd and RR at z%d done"%ind
        for ps in [0, 15, 25, 35, 45]:
            for j in linds:
                halos = lpath%(ps, ind, ps, j, ind)
                R, xihm = calc_hmcf(halos, RR, DdRd, Dd, Rh, hmcf_savepath%(ps,i,j))
                print "HMCF ps%d z%d, l%d complete"%(ps,ind,j)
                continue #end j
            continue #end ps
        continue #end i,ind
    print "HMCFs created"

    """
    for ps in [15, 25, 35]:
        mean_masses = np.genfromtxt("L_ps%d_masses.txt"%ps)
        for i,ind in zip(range(len(inds)), inds):
            z = zs[i]
            for j in linds:
                Mass = mean_masses[i,j]
                R, xihm = np.loadtxt(hmcf_savepath%(ps,i,j))
                hminds = np.invert(np.isnan(xihm)+(xihm==1.0))
                R = R[hminds]
                xihm = xihm[hminds]
                R, DS = calc_DS(R, xihm, Mass, z, ds_savepath%(ps,i,j))
                print "DeltaSigma ps%d z%d, l%d complete"%(ps,ind,j)
                continue #end j
            continue #end i,ind
        continue #end i,ind
    print "DeltaSigmas created"
    """

    """
    for ps in [15, 25, 35]:
        for i,ind in zip(range(len(inds)), inds):
            z = zs[i]
            for j in linds:
                R, DS = np.loadtxt(ds_savepath%(ps,i,j)).T
                cov = np.genfromtxt(covbase%(j, covinds[i]))
                DSsave  = DSdatabase%(ps,i,j)
                covsave = covdatabase%(ps,i,j)
                create_data_vector(R, DS, cov, z, DSsave, covsave)
                print "Data vector created for ps%d z%d, l%d"%(ps, ind, j)
                continue  #end j
            continue #end i,ind
        continue #end ps
    print "Data vectors created"
    """
