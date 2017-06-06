"""
This is a file that combines all the small scripts I used earlier.
This is being written so that I can use general paths to the catalogs of 
interest isntead of hard-coding in paths to certain catalogs.
"""
import os, sys
import numpy as np
import pygadgetreader as pgr
from Corrfunc.theory.xi import xi
from Corrfunc.theory.DD import DD
from Corrfunc.utils import convert_3d_counts_to_cf

inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
halobase="/calvin1/tmcclintock/fox_data/richness_halos/rich_snapdir_ps%d_%03d/"
mpath = halobase+"/mass_halos_ps%d_m%d_%03d.txt" #Path to mass split halos
lpath = halobase+"/richness_halos_ps%d_l%d_%03d.txt" #Path to lam split halos
lam_edges = [5, 10, 14, 20, 30, 45, 60, np.inf]
lM_edges = [13.1, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.5, 15.0]#, 16.0]

dmpath = "/calvin1/tmcclintock/down_sampled_snapshots/snapdir_%03d/snapshot_%03d_z%s_down1000"

boxsize = 1050.0
nthreads = 8
nbins = 50
bins = np.logspace(np.log10(0.01), np.log10(150.0), nbins+1)

def calc_hhcf(pscatter, domass = False):
    print "Calculating hhcf for pscatter=%d"%pscatter
    os.system("mkdir -p txt_files/xihh_ps%d"%pscatter)
    for i, red in zip(inds, zs):
        for j in range(len(lam_edges)-1):
            X, Y, Z, Np, M, Rich = np.genfromtxt(lpath%(pscatter, i, pscatter, j,i), unpack=True)
            result = xi(boxsize, nthreads, bins, X, Y, Z)
            np.savetxt("txt_files/xihh_ps%d/xihh_ps%d_z%.2f_l%d.txt"%(pscatter, pscatter, red, j), result)
            print "xihh for z%.2f l%d"%(red, j)
            continue
        continue
    return

def calc_hmcf(pscatter, domass = False):
    print "Calculating hmcf for pscatter=%d"%pscatter
    os.system("mkdir -p txt_files/xihm_ps%d"%pscatter)
    for i, red, zstring in zip(inds, zs, zstrings):
        Xd, Yd, Zd = pgr.readsnap(dmpath%(i, i, zstring), 'pos', 'dm').T
        Xd = Xd.astype('float64')
        Yd = Yd.astype('float64')
        Zd = Zd.astype('float64')
        N_dm = len(Xd)
        N_rand = int(N_dm*1.5)
        Xr1 = np.random.uniform(0, boxsize, N_rand)
        Yr1 = np.random.uniform(0, boxsize, N_rand)
        Zr1 = np.random.uniform(0, boxsize, N_rand)
        Xr2 = np.random.uniform(0, boxsize, N_rand)
        Yr2 = np.random.uniform(0, boxsize, N_rand)
        Zr2 = np.random.uniform(0, boxsize, N_rand)
        # do the RR and the DdRd term
        RhRd = DD(0, nthreads, bins, Xr1, Yr1, Zr1, X2=Xr2, Y2=Yr2, Z2=Zr2)
        DdRd = DD(0, nthreads, bins, Xd, Yd, Zd, X2=Xr2, Y2=Yr2, Z2=Zr2)
        for j in range(len(lam_edges)-1):
            Xh, Yh, Zh, Np, M, Rich = np.genfromtxt(lpath%(pscatter, i, pscatter, j,i), unpack=True)
            DhDd = DD(0, nthreads, bins, X1=Xh, Y1=Yh, Z1=Zh, 
                      X2=Xd, Y2=Yd, Z2=Zd)
            DhRh = DD(0, nthreads, bins, Xh, Yh, Zh,
                      X2=Xr1, Y2=Yr1, Z2=Zr1)
            hmcf = convert_3d_counts_to_cf(N_h, N_dm, N_rand, N_rand,
                                           DhDd, DhRh, DdRd, RhRd)
            np.savetxt("txt_files/xihm_ps%d/xihm_ps%d_z%.2f_l%d.txt"%(pscatter, pscatter, red, j), hmcf)
            print "xihm for z%.2f l%d"%(red, j)
            continue
        continue
    return

if __name__ == "__main__":
    pscatter = 25
    #calc_hhcf(pscatter)
    calc_hmcf(pscatter)
