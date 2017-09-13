import os, sys
sys.path.insert(0,"./src/")
from likelihoods import *
from get_data import *
import matplotlib.pyplot as plt

#Fox cosmology
h = 0.670435

def plot_data(ps, zi, lj, use_y1=True):
    z = zs[zi]
    klin, Plin, knl, Pmm = get_P(zi)
    cosmo, input_params = get_cosmo_and_params()
    input_params["R_bin_min"] = 0.0323*(h*(1+z)) #Mpc/h comoving
    input_params["R_bin_max"] = 30.0*(h*(1+z)) #Mpc/h comoving
    extras = [klin, knl, Plin, Pmm, cosmo, input_params]
    R, DS, cov, cut = get_data_and_cov(ps, zi, lj, use_y1, nocut=True)
    DSe = np.sqrt(np.diag(cov))
    bad = R < 0.2 #Mpc
    gud = R > 0.2 #Mpc
    plt.errorbar(R[bad], DS[bad], DSe[bad], ls='', marker='.', mfc='w', c='b')
    plt.errorbar(R[gud], DS[gud], DSe[gud], ls='', marker='.', c='b')
    plt.xscale('log')
    plt.yscale('log')
    #plt.show()

def plot_bf(ps, zi, li, use_y1=True):
    print "working on this"

if __name__ == "__main__":
    plot_data(0, 0, 0)
    plot_bf(0, 0, 0)
    plt.show()
