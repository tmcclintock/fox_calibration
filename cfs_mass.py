"""
Run corrfunc on the fox halos that are split by mass.

Richness bin edges: 5, 10, 14, 20, 30, 45, 60, inf for richness.
They are given by mass below by lM_edges.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
import pygadgetreader as pgr
from Corrfunc.theory.xi import xi
from Corrfunc.theory.DD import DD
from Corrfunc.utils import convert_3d_counts_to_cf

halopath = "/calvin1/tmcclintock/fox_data/richness_halos/mass_snapdir_%03d/mass_halos_m%d_%03d.txt"
dmpath = "/calvin1/tmcclintock/down_sampled_snapshots/snapdir_%03d/snapshot_%03d_z%s_down100"

boxsize = 1050.0
nthreads = 8
nbins = 50
bins = np.logspace(np.log10(0.01), np.log10(150.0), nbins+1)
R = (bins[:-1]+bins[1:])/2.

dohh = False
calchh = False
dohm = True
calchm = True

inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
lM_edges = [13.1, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.5, 15.0]#, 16.0]
minds = range(0, len(lM_edges)-1)
c = np.linspace(1.0, 0.3, len(minds))

cmaps = ['Reds','Oranges','Greens', 'Blues']

if dohh:
    for i in range(len(inds)):
        index = inds[i]
        z = zs[i]
        for j in minds:
            if calchh:
                print "finding HHCF at z=%.2f m=%d"%(z,j)
                X, Y, Z, Np, M, Rich = np.genfromtxt(halopath%(index,j,index), unpack=True)
                result = xi(boxsize, nthreads, bins, X, Y, Z)
                np.savetxt("txt_files/mass_txt_files/hhcf_z%.2f_m%d.txt"%(z,j), result)
        for j in minds:
            result = np.loadtxt("txt_files/mass_txt_files/hhcf_z%.2f_m%d.txt"%(z,j)).T
            hhcf = result[3]
            cmap = plt.get_cmap(cmaps[i])
            plt.loglog(R, hhcf, c=cmap(c[j]), label="z=%.2f m%d"%(z,j))
        plt.legend(loc=0)
        plt.xlabel(r"$R\ [{\rm Mpc/h}]$", fontsize=24)
        plt.ylabel(r"$\xi_{\rm hh}$", fontsize=24)
        plt.subplots_adjust(bottom=0.15)
        plt.gcf().savefig("figures/hhcf_masses_z%0.2f.png"%z)
        plt.show()
        plt.clf()

if dohm:
    for i in range(len(inds)):
        #if i == 0: continue
        cmap = plt.get_cmap(cmaps[i])
        index = inds[i]
        z = zs[i]
        zstring = zstrings[i]
        if calchm:
            print "finding HMCF at z=%.2f"%z
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
            print "starting on RR and DdRd"
            RhRd = DD(0, nthreads, bins, Xr1, Yr1, Zr1,
                      X2=Xr2, Y2=Yr2, Z2=Zr2)
            print "starting DdRd"
            DdRd = DD(0, nthreads, bins, Xd, Yd, Zd,
                      X2=Xr2, Y2=Yr2, Z2=Zr2)
            print "finished with RR and DdRd"

            for j in minds:
                print "at m%d"%j
                Xh, Yh, Zh, Np, M, Rich = np.genfromtxt(halopath%(index,j,index), unpack=True)
                N_h = len(Xh)
                DhDd = DD(0, nthreads, bins, X1=Xh, Y1=Yh, Z1=Zh, 
                          X2=Xd, Y2=Yd, Z2=Zd)
                DhRh = DD(0, nthreads, bins, Xh, Yh, Zh,
                          X2=Xr1, Y2=Yr1, Z2=Zr1)
                hmcf = convert_3d_counts_to_cf(N_h, N_dm, N_rand, N_rand,
                                               DhDd, DhRh, DdRd, RhRd)
                np.savetxt("txt_files/mass_txt_files/hmcf_z%.2f_m%d.txt"%(z,j),hmcf)
        for j in minds:
            hmcf = np.loadtxt("txt_files/mass_txt_files/hmcf_z%.2f_m%d.txt"%(z,j))
            plt.loglog(R, hmcf, label="z=%.2f m%d"%(z,j))
        plt.legend(loc=0)
        plt.xlabel(r"$R\ [{\rm Mpc/h}]$", fontsize=24)
        plt.ylabel(r"$\xi_{\rm hm}$", fontsize=24)
        plt.subplots_adjust(bottom=0.15)
        plt.gcf().savefig("figures/hmcf_masses_z%0.2f.png"%z)
        #plt.show()
        plt.clf()
