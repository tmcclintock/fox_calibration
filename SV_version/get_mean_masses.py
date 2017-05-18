"""
Calculate the mean masses.
"""

import numpy as np

#Dimensions of the output array
nz = 4
nl = 7
inds = [6,7,8,9]

mass = np.zeros((nz, nl))

inpath = ""
halopath = "/calvin1/tmcclintock/fox_data/richness_halos/rich_snapdir_%03d/richness_halos_l%d_%03d.txt"


for i in range(len(inds)):
    index = inds[i]
    for j in range(nl):
        X, Y, Z, Np, M, Rich = np.genfromtxt(halopath%(index,j,index), 
                                             unpack=True)
        mass[i,j] = np.mean(M)
        print "Have mass for %d %d"%(index,j)

print mass
np.savetxt("txt_files/mean_masses.txt", mass)
