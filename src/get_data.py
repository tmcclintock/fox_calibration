"""
This contains routines to get the correct data out.
"""

import numpy as np

#Various paths to the data files
inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
DSdatabase      = "/calvin1/tmcclintock/DES_Y1_data/calibration_data/dss/ds_ps%d_z%d_l%d.txt"
covdatabase     = "/calvin1/tmcclintock/DES_Y1_data/calibration_data/covs/cov_ps%d_z%d_l%d.txt"
svcovdatabase   = "/calvin1/tmcclintock/DES_Y1_data/calibration_data/svcovs/svcov_ps%d_z%d_l%d.txt"

def get_data_and_cov(ps, i, j, use_y1=True, nocut=False):
    DSpath  = DSdatabase%(ps, i, j)
    if use_y1: covpath = covdatabase%(ps, i, j)
    else: covpath = svcovdatabase%(ps, i, j)
    R, DS = np.loadtxt(DSpath).T
    cov = np.loadtxt(covpath)
    if nocut: cut = R > 0.0 #everything
    else: cut = R>0.2 #Mpc
    cov = cov[cut]
    cov = cov[:,cut]
    DS = DS[cut]
    R = R[cut]
    return R, DS, cov, cut

#Input parameters for the DeltaSigma code
#Fox cosmology
h = 0.670435
cosmo = {"h":h,"om":0.31834,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]
input_params = {"NR":300,"Rmin":0.01,
                "Rmax":200.0,"Nbins":15,"delta":200,
                "Rmis":0.0, "fmis":0.0,
                "miscentering":0,"averaging":1,"single_miscentering":0}
#Dictionary for colossus
colcos = {"H0":cosmo['h']*100.,"Om0":cosmo['om'], 
          'Ob0': 0.049017, 'sigma8': 0.83495, 'ns': 0.96191, 'flat':True}


def get_cosmo_and_params():
    return cosmo, input_params
def get_colcos():
    return colcos

#The power spectrum
klin = np.loadtxt("txt_files/P_files/k.txt")
def get_P(i, use_old_P=False):
    z = zs[i]
    #use_old_P can let us use the non-Takahashi P_mm
    Plin = np.loadtxt("txt_files/P_files/Plin_z%.2f.txt"%z)
    if use_old_P: knl  = np.loadtxt("txt_files/P_files/knl.txt")
    else: knl = np.loadtxt("txt_files/P_files/k.txt")
    if use_old_P: Pmm  = np.loadtxt("txt_files/P_files/Pnl_old_z%.2f.txt"%z)
    else: Pmm  = np.loadtxt("txt_files/P_files/Pnl_z%.2f.txt"%z)
    return klin, Plin, knl, Pmm

#True masses
def get_true_masses(ps, M_split=False):
    return np.loadtxt("L_ps%d_masses.txt"%ps)

#best fit results
def get_bf_results(ps, use_y1=True):
    if use_y1: bfpath = "output_files/mass_fits/bf_y1_masses_ps%d.txt"%ps
    else: bfpath = "output_files/mass_fits/bf_sv_masses_ps%d.txt"%ps
    return np.loadtxt(bfpath)

