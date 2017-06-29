"""
This contains the post-processing steps
"""
import numpy as np
import matplotlib.pyplot as plt

bfpath  = "output_files/mass_fits/bf_cal_ps%d.txt"
calpath = "output_files/mass_fits/mcmc_cal_ps%d.txt"
errpath = "output_files/mass_fits/mcmc_calerr_ps%d.txt"
richpath = "L_ps%d_richnesses.txt"

zs = [1.0, 0.5, 0.25, 0.0]
fig, axarr = plt.subplots(len(zs), sharex=True, sharey=True)
if __name__ == "__main__":
    pss = [15, 25, 35]

    for ps,j in zip(pss,range(len(pss))):
        bf  = np.loadtxt(bfpath%ps)
        cal = np.loadtxt(calpath%ps)
        err = np.loadtxt(errpath%ps)
        lams = np.genfromtxt(richpath%ps)
        for i in range(len(zs)):
            hi = lams[i]>20
            lo = lams[i]<20
            z = zs[i]
            axarr[i].plot(lams[i], bf[i])
            axarr[i].errorbar(lams[i,hi], cal[i,hi], err[i,hi], marker='o',ls='',c='C%d'%j, alpha=0.5)
            axarr[i].errorbar(lams[i,lo], cal[i,lo], err[i,lo], marker='o',mfc='white',ls='',c='C%d'%j, alpha=0.5)
            if j==0:
                axarr[i].text(60, 1.15, r"$z=%.2f$"%zs[i], fontsize=16)

    plt.show()
