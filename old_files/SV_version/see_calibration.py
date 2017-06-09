"""
See the calibration curves.

Note: we only have best fit masses so far.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
plt.rc("text", usetex=True, fontsize=24)
plt.rc("errorbar", capsize=3)

use_BF = False
if use_BF:
    bfs  = np.loadtxt("txt_files/SV_BF_masses.txt")
else:
    bfs  = np.loadtxt("txt_files/SV_mcmc_masses.txt")
    errs = np.loadtxt("txt_files/SV_mcmc_errs.txt")
mass = np.loadtxt("txt_files/SV_masses.txt")
C = mass/10**bfs
if not use_BF:
    Cerr = np.log(10) * errs * C

zs = [1.0, 0.5, 0.25, 0.0]

fig, axarr = plt.subplots(4, sharex=True, sharey=True)
for i in range(len(zs)):
    z = zs[i]
    print C[i]
    if use_BF:
        axarr[i].plot(bfs[i], C[i], marker='o', ls='')
    else:
        axarr[i].errorbar(bfs[i], C[i], Cerr[i], marker='o', ls='')
    axarr[i].axhline(y=1.05, c='k', ls='--', lw=1)

for i in range(len(zs)):
    ylim = .32#max(axarr[i].get_ylim())-1 #% above 1
    #axarr[i].set_ylim(1-ylim, 1+ylim)
    axarr[i].set_ylim(1, 1.3)
    axarr[i].tick_params(labelsize=16)
    axarr[i].text(14, 1.15, r"$z=%.2f$"%zs[i], fontsize=16)
    #nbins = len(axarr[i].get_yticklabels())
    #axarr[i].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

plt.subplots_adjust(hspace=0.1, bottom=0.2, left=0.2)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel(r"$M_{\rm obs}\ [{\rm M_\odot/h}]$")
plt.ylabel(r"$\cal{C}=\frac{M_{\rm true}}{M_{\rm obs}}$")
plt.show()
