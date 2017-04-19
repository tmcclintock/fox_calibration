"""
Run corrfunc on the fox halos to produce correlation functions.

First I calculate halo-halo correlation functions for
all 
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
import pygadgetreader as pgr
from Corrfunc.theory.xi import xi
from Corrfunc.theory.DD import DD
from Corrfunc.utils import convert_3d_counts_to_cf

halopath = "/calvin1/tmcclintock/fox_data/richness_halos/rich_snapdir_%03d/reduced_richness_halos_%03d"
dmpath = "/calvin1/tmcclintock/down_sampled_snapshots/snapdir_%03d/snapshot_%03d_z%s_down100"

boxsize = 1050.0
nthreads = 8
nbins = 100
bins = np.logspace(np.log10(0.01), np.log10(150.0), nbins+1)
R = (bins[:-1]+bins[1:])/2.

dohh = False
calchh = False
dohm = True
calchm = False

inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]

if dohh:
    if calchh:
        for i in range(len(inds)):
            index = inds[i]
            z = zs[i]
            print "finding HHCF at z=%.2f"%z
            X, Y, Z, Np, M, Rich = np.genfromtxt(halopath%(index,index), unpack=True)
            print len(X)
            result = xi(boxsize, nthreads, bins, X, Y, Z)
            np.savetxt("txt_files/hhcf_z%.2f.txt"%z, result)
    for i in range(len(inds)):
        index = inds[i]
        z = zs[i]
        result = np.loadtxt("txt_files/hhcf_z%.2f.txt"%z).T
        r = np.mean(result[:2], 0)
        xi = result[3]
        plt.loglog(r, xi, label="z=%.2f"%z)
    plt.legend(loc=0)
    plt.xlabel(r"$R\ [{\rm Mpc/h}]$", fontsize=24)
    plt.ylabel(r"$\xi_{\rm hh}$", fontsize=24)
    plt.subplots_adjust(bottom=0.15)
    plt.show()
    plt.clf()

if dohm:
    if calchm:
        for i in range(len(inds)):
            index = inds[i]
            z = zs[i]
            if z > 0.0: continue
            zstring = zstrings[i]
            print "finding HMCF at z=%.2f"%z
            Xh, Yh, Zh, Np, M, Rich = np.genfromtxt(halopath%(index,index), unpack=True)
            N_h = len(Xh)
            Xd, Yd, Zd = pgr.readsnap(dmpath%(index,index,zstring), 'pos', 'dm').T
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

            DhDd = DD(0, nthreads, bins, X1=Xh, Y1=Yh, Z1=Zh, 
                      X2=Xd, Y2=Yd, Z2=Zd)
            DhRh = DD(0, nthreads, bins, Xh, Yh, Zh,
                      X2=Xr1, Y2=Yr1, Z2=Zr1)
            DdRd = DD(0, nthreads, bins, Xd, Yd, Zd,
                      X2=Xr2, Y2=Yr2, Z2=Zr2)
            RhRd = DD(0, nthreads, bins, Xr1, Yr1, Zr1,
                      X2=Xr2, Y2=Yr2, Z2=Zr2)
            hmcf = convert_3d_counts_to_cf(N_h, N_dm, N_rand, N_rand,
                                           DhDd, DhRh, DdRd, RhRd)
            np.savetxt("txt_files/hmcf_z%.2f.txt"%z,hmcf)
    for i in range(len(inds)):
        index = inds[i]
        z = zs[i]
        xi = np.loadtxt("txt_files/hmcf_z%.2f.txt"%z)
        print R.shape, xi.shape
        plt.loglog(R, xi, label="z=%.2f"%z)
    plt.legend(loc=0)
    plt.xlabel(r"$R\ [{\rm Mpc/h}]$", fontsize=24)
    plt.ylabel(r"$\xi_{\rm hm}$", fontsize=24)
    plt.subplots_adjust(bottom=0.15)
    plt.show()
    plt.clf()
