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
mass = np.loadtxt("txt_files/true_masses.txt")
if use_BF:
    bfs  = np.loadtxt("txt_files/BF_masses.txt")
else:
    bfs  = np.loadtxt("txt_files/mcmc_masses.txt")
    errs = np.loadtxt("txt_files/mcmc_errs.txt")
C = mass/10**bfs
if not use_BF:
    Cerr = np.log(10) * errs * C

lams = np.loadtxt("txt_files/true_richnesses.txt")
zs = [1.0, 0.5, 0.25, 0.0]

def calc_cal(lam, z):
    C0, a , b =  [1.06569352, -0.00416351, 0.01785429]
    dC0, da, db = [ 0.02810846, 0.08770361, 0.05685181]
    Cbf = C0*((1+z)/1.5)**a*(lam/30.)**b
    Cbfe = Cbf**2/C0**2*dC0**2 + Cbf**2*np.log(lam/30.)**2*(da**2 + db**2)
    return Cbf, Cbfe

def get_cal_curve(z):
    """
    Return four arrays:
    - richness values
    - C best fit values
    - +/- error bars
    """
    lam = np.linspace(4, 100, 100)
    Cbf, Cbfe = calc_cal(lam, z)
    return [lam, Cbf, Cbfe]

fig, axarr = plt.subplots(4, sharex=True, sharey=True)
for i in range(len(zs)):
    z = zs[i]
    hi = lams[i]>20.
    lo = lams[i]<20.
    if use_BF:
        axarr[i].plot(lams[i,hi], C[i,hi], marker='o', ls='')
        axarr[i].plot(lams[i,lo], C[i,lo], marker='o', mfc='white', ls='')
    else:
        axarr[i].errorbar(lams[i,hi], C[i,hi], Cerr[i,hi], marker='o', ls='', c='C0')
        axarr[i].errorbar(lams[i,lo], C[i,lo], Cerr[i,lo], marker='o', mfc='white', ls='',c='C0')

    print calc_cal(lams[i,hi], z)[0]
    lamx, Cbf, Cbfe = get_cal_curve(z)
    axarr[i].fill_between(lamx, Cbf-Cbfe, Cbf+Cbfe, color="r",alpha=0.2, zorder=-1)
    #axarr[i].axhline(y=1.05, c='k', ls='--', lw=1)

for i in range(len(zs)):
    ylim = .32#max(axarr[i].get_ylim())-1 #% above 1
    #axarr[i].set_ylim(1-ylim, 1+ylim)
    axarr[i].set_ylim(1, 1.3)
    axarr[i].set_xlim(4, 100)
    axarr[i].tick_params(labelsize=16)
    axarr[i].text(60, 1.15, r"$z=%.2f$"%zs[i], fontsize=16)
    #nbins = len(axarr[i].get_yticklabels())
    #axarr[i].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

plt.subplots_adjust(hspace=0.1, bottom=0.2, left=0.2)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel(r"$\lambda$")
plt.ylabel(r"${\cal C}=\frac{M_{\rm true}}{M_{\rm obs}}$")
plt.show()
