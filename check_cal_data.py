import os, sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "./src/")
import clusterwl
from CF_functions import *
from DS_functions import *
import helper_functions as HF
plt.rc("text", usetex=True)
plt.rc("font", size=18)
cosmo = HF.get_cosmo()
h = cosmo['h']
om = cosmo['om']
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

zs = [1.0, 0.5, 0.25, 0.0]
hmcf_savepath = "output_files/hmcf/hmcf_ps%d_z%d_l%d.txt"
ds_savepath   = "output_files/ds/ds_ps%d_z%d_l%d.txt"
svDSdatabase = "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/svds_ps%d_z%d_l%d.txt"
y1DSdatabase = "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/y1ds_ps%d_z%d_l%d.txt"
svMbfpath = "output_files/mass_fits/bf_sv_masses_ps%d.txt"
y1Mbfpath = "output_files/mass_fits/bf_y1_masses_ps%d.txt"
svCbfpath = "output_files/mass_fits/bf_sv_conc_ps%d.txt"
y1Cbfpath = "output_files/mass_fits/bf_y1_conc_ps%d.txt"


def get_xi_model(M, c, z):
    k, Plin, Pnl = HF.get_P(z)
    Rmodel = np.logspace(-2, 3, num=1000, base=10) 
    xi_mm = clusterwl.xi.xi_mm_at_R(Rmodel, k, Pnl)
    xi_nfw   = clusterwl.xi.xi_nfw_at_R(Rmodel, M, c, om)
    bias = clusterwl.bias.bias_at_M(M, k, Plin, om)
    xi_2halo = clusterwl.xi.xi_2halo(bias, xi_mm)
    xi_hm    = clusterwl.xi.xi_hm(xi_nfw, xi_2halo)
    return Rmodel, xi_nfw, xi_2halo, xi_hm

def get_ds_model(M, c, z):
    Rmodel, xi_nfw, xi_2halo, xi_hm = get_xi_model(M, c, z)
    Rp = np.logspace(-2, 2.4, 1000, base=10) #Mpc/h
    Sigma  = clusterwl.deltasigma.Sigma_at_R(Rp, Rmodel, xi_hm, M, c, om)
    DeltaSigma = clusterwl.deltasigma.DeltaSigma_at_R(Rp, Rp, Sigma, M, c, om)
    return Rp, DeltaSigma

def check_hmcfs(ps, i, j):
    R, xihm = np.loadtxt(hmcf_savepath%(ps,i,j)) #first the data
    z = zs[i]
    mean_masses = np.genfromtxt("txt_files/L_ps%d_masses.txt"%ps)
    lM = np.log10(mean_masses[i, j])
    M = 10**lM
    c = get_conc(M, z)
    Rmodel, xi_nfw, xi_2halo, xi_hm_model = get_xi_model(M, c, z)
    xispl = IUS(Rmodel, xi_hm_model)
    pdiff = (xihm - xispl(R))/xispl(R)

    f, axarr = plt.subplots(2, sharex=True)
    
    axarr[0].loglog(Rmodel, xi_nfw, c='red', ls='--', label=r'$1-{\rm halo}')
    axarr[0].loglog(Rmodel, xi_2halo, c='orange', ls='--', label=r'$2-{\rm halo}$')
    axarr[0].loglog(R, np.fabs(xihm), c='blue', label=r'${\rm Sim}$')
    axarr[1].plot(R, pdiff, c='k')
    axarr[1].axhline(0, ls='--', c='k')

    leg = axarr[0].legend(loc=0, fancybox=True, frameon=False, fontsize=16)
    plt.xlim(1e-1, 4e1)
    axarr[0].set_ylim(1e-2, 2e4)
    P = 0.30
    #axarr[1].set_ylim(-0.33, 0.33)
    axarr[1].set_ylim(-P, P)
    plt.xlabel(r"$R\ [{\rm Mpc}/h]$")
    axarr[0].set_ylabel(r"${\rm Halo-matter\ clustering}$")
    axarr[1].set_ylabel(r"$\%\ {\rm Diff}$")
    plt.subplots_adjust(bottom=0.15, left=0.18, hspace=0)
    #plt.gcf().savefig("xihm.eps")
    plt.gcf().savefig("xihm.pdf")
    plt.show()

def check_ds(ps, i, j, bf=False, y1=True):
    z = zs[i]
    mean_masses = np.genfromtxt("txt_files/L_ps%d_masses.txt"%ps)
    lM = np.log10(mean_masses[i, j])
    M = 10**lM
    c = get_conc(M, z)
    Rp, DeltaSigma = get_ds_model(M, c, z)

    if bf:
        if y1: 
            Mbf = np.loadtxt(y1Mbfpath%ps)[i, j]
            cbf = np.loadtxt(y1Cbfpath%ps)[i, j]
        else:  
            Mbf = np.loadtxt(svMbfpath%ps)[i, j]
            cbf = np.loadtxt(svCbfpath%ps)[i, j]
        Rbf, DSbf = get_ds_model(Mbf, cbf, z)
        plt.loglog(Rbf, DSbf, c='k', label=r"$\Delta\Sigma(M_{best})$")

    Rpd, DS = np.loadtxt(ds_savepath%(ps,i,j), unpack=True)
    #These data are in Msun, Mpc physical
    Rmidsv, svDS = np.loadtxt(svDSdatabase%(ps, i, j), unpack=True)
    Rmidy1, y1DS = np.loadtxt(y1DSdatabase%(ps, i, j), unpack=True)

    plt.loglog(Rp, DeltaSigma, c='r', ls='--', label=r"$\Delta\Sigma(M_{true})$")
    plt.loglog(Rpd, DS, label='sim')
    
    #Convert the data to Msun/h comoving
    #plt.loglog(Rmidsv*h*(1+z), svDS/(h*(1+z)**2), label='sv sim', marker='.', ls='')
    #plt.loglog(Rmidy1*h*(1+z), y1DS/(h*(1+z)**2), label='y1 sim', marker='.', ls='')
    plt.legend(loc = 0)
    plt.show()


if __name__ == "__main__":
    ps, zi, lj = 0, 3, 4
    check_hmcfs(ps, zi, lj)
    #check_ds(ps, zi, lj, bf=False, y1=False)
