import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)

nbins = 50
bins = np.logspace(np.log10(0.01), np.log10(150.0), nbins+1)
R = (bins[:-1]+bins[1:])/2.

inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
c = np.linspace(1.0, 0.0, len(linds))
cmaps = ['Reds','Oranges','Greens', 'Blues']

for i in range(len(inds)):
    if i >=2: continue
    cmap = plt.get_cmap(cmaps[i])
    index = inds[i]
    z = zs[i]
    zstring = zstrings[i]
    for j in linds:
        hmcf = np.loadtxt("txt_files/richness_txt_files/hmcf_z%.2f_l%d.txt"%(z,j))
        plt.loglog(R, hmcf, label="z=%.2f l%d"%(z,j))
    plt.legend(loc=0)
    plt.xlabel(r"$R\ [{\rm Mpc/h}]$", fontsize=24)
    plt.ylabel(r"$\xi_{\rm hm}$", fontsize=24)
    plt.subplots_adjust(bottom=0.15)
    plt.gcf().savefig("figures/hmcf_richnesses_z%0.2f.png"%z)
    plt.show()
    plt.clf()
