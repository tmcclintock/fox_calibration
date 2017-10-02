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
import clusterwl
from CF_functions import *
from DS_functions import *
import helper_functions as HF

cosmo = HF.get_cosmo()
h = cosmo['h']

inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
halobase= "/calvin1/tmcclintock/fox_data/richness_halos/rich_snapdir_ps%d_%03d/"
mpath = halobase+"/mass_halos_ps%d_m%d_%03d.txt" #Path to mass split halos
lpath = halobase+"/richness_halos_ps%d_l%d_%03d.txt" #Path to lam split halos
lam_edges = [5, 10, 14, 20, 30, 45, 60, np.inf]
lM_edges = [13.1, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.5, 15.0]#, 16.0]

dmpath_base = "/calvin1/tmcclintock/down_sampled_snapshots/snapdir_%03d/snapshot_%03d_z%s_down100"
hhcf_savepath = "output_files/hhcf/hhcf_ps%d_z%d_l%d.txt"
hmcf_savepath = "output_files/hmcf/hmcf_ps%d_z%d_l%d.txt"
ds_savepath   = "output_files/ds/ds_ps%d_z%d_l%d.txt"

def make_HHCFs():
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

def make_HMCFs():
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

fory1 = True
if fory1: DSdatabase = "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/y1ds_ps%d_z%d_l%d.txt"
else:     DSdatabase = "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/svds_ps%d_z%d_l%d.txt"

def make_DSs(ps):
    mean_masses = np.genfromtxt("txt_files/L_ps%d_masses.txt"%ps)
    for i,ind in zip(range(len(inds)), inds):
        z = zs[i]
        k, Plin, Pnl = HF.get_P(z)
        Rmodel = np.logspace(-2, 3, num=1000, base=10) 
        xi_mm = clusterwl.xi.xi_mm_at_R(Rmodel, k, Pnl)

        for j in linds:
            Mass = mean_masses[i,j]
            R, xihm = np.loadtxt(hmcf_savepath%(ps,i,j))
            hminds = np.invert(np.isnan(xihm)+(xihm==1.0))
            R = R[hminds]
            xihm = xihm[hminds]
            Rfull, DS, Rmid, aDS = calc_DS(R, xihm, Mass, z, Rmodel, xi_mm, k, Plin, Pnl, dssave=ds_savepath%(ps, i, j), avsave = DSdatabase%(ps, i, j), y1_binning=fory1)
            print "DeltaSigma ps%d z%.2f, l%d complete"%(ps,z,j)
            continue #end j
    print "DeltaSigmas created for ps%d"%ps

y1covbase = "/calvin1/tmcclintock/DES_DATA_FILES/y1_data_files/blinded_tamas_files/full-mcal-raw_y1subtr_l%d_z%d_dst_cov.dat"
svcovbase = "/calvin1/tmcclintock/DES_DATA_FILES/sv_data_files/cov_t_z%d_l%d.dat"

y1covoutbase= "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/y1cov_z%d_l%d.txt"
svcovoutbase= "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/svcov_z%d_l%d.txt"
covzinds = [2, 2, 2, 1] #z indices for covariance matrices to use for each snap
svcovzinds  = [2, 1, 0, 0] #z indices for matchin SV covariance matrices to each snap
svcovlinds = [0, 1, 2, 3, 4, 4, 4] #l indices for matching SV covariance matrices to each snap

def save_cov():
    for i,ind in zip(range(len(inds)), inds):
        for j in linds:
            y1cov   = np.genfromtxt(y1covbase%(j, covzinds[i]))
            svcov = np.genfromtxt(svcovbase%(svcovzinds[i], svcovlinds[j]))
            y1covsave   = y1covoutbase%(i,j)
            svcovsave = svcovoutbase%(i,j)
            np.savetxt(y1covsave, y1cov)
            np.savetxt(svcovsave, svcov)
            print "Covs saved for z%d, l%d"%(ind, j)
            continue  #end j
        continue #end i,ind
    print "Covariances saved"

if __name__ == "__main__":
    make_DSs(0)
    make_DSs(15)
    make_DSs(25)
    make_DSs(35)
    #save_cov()
