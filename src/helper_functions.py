"""
This contains routines to get the correct data out.
"""

import numpy as np

#Various paths to the data files
inds = [6,7,8,9]
zs = np.array([1.0, 0.5, 0.25, 0.0])
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
y1DSdatabase      = "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/y1ds_ps%d_z%d_l%d.txt"
svDSdatabase      = "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/svds_ps%d_z%d_l%d.txt"
y1covdatabase     = "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/y1cov_z%d_l%d.txt"
svcovdatabase   = "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/svcov_z%d_l%d.txt"

def get_zs():
    return zs, zstrings

def get_lams(ps):
    return np.loadtxt("txt_files/L_ps%d_richnesses.txt"%ps)

def get_true_M(ps):
    return np.loadtxt("txt_files/L_ps%d_masses.txt"%ps)

def get_Redges(use_y1=True):
    #The bin edges in Mpc physical
    Nbins = 15
    if use_y1: return np.logspace(np.log10(0.0323), np.log10(30.), num=Nbins+1)
    else: return np.logspace(np.log10(0.02), np.log10(30.), num=Nbins+1) #use_sv

def get_data_and_cov(ps, i, j, use_y1=True, nocut=False):
    if use_y1: 
        DSpath  = y1DSdatabase%(ps, i, j)
        covpath = y1covdatabase%(i, j)
    else: 
        DSpath  = svDSdatabase%(ps, i, j)
        covpath = svcovdatabase%(i, j)
    R, DS = np.loadtxt(DSpath).T
    cov = np.loadtxt(covpath)
    if nocut: return R, DS, np.linalg.inv(cov), cov, inds
    cuts = [0.2, 999]
    if use_y1==False and i > 1: cuts[1] = 21.5
    inds = (R > cuts[0])*(R < cuts[1]) #Mpc
    cov = cov[inds]
    cov = cov[:,inds]
    DS = DS[inds]
    R = R[inds]
    return R, DS, np.linalg.inv(cov), cov, inds

#Input parameters for the DeltaSigma code
#Fox cosmology
h = 0.670435
cosmo = {"h":h,"om":0.31834,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]
#Dictionary for colossus
colcos = {"H0":cosmo['h']*100.,"Om0":cosmo['om'], 
          'Ob0': 0.049017, 'sigma8': 0.83495, 'ns': 0.96191, 'flat':True}

def get_cosmo():
    return cosmo
def get_cosmo_and_params():
    return cosmo, input_params
def get_colcos():
    return colcos

#The power spectrum
P_path = "/calvin1/tmcclintock/DES_DATA_FILES/fox_data_files/P_files/"
def get_P(z):
    k = np.loadtxt(P_path + "k.txt")
    Plin = np.loadtxt(P_path + "plin_z%.2f.txt"%z)
    Pnl  = np.loadtxt(P_path + "pnl_z%.2f.txt"%z)
    return k, Plin, Pnl

#best fit results
def get_bf_results(ps, use_y1=True):
    if use_y1: bfpath = "output_files/mass_fits/bf_y1_masses_ps%d.txt"%ps
    else: bfpath = "output_files/mass_fits/bf_sv_masses_ps%d.txt"%ps
    return np.loadtxt(bfpath)

def get_model_defaults(h):
    #Dictionary of default starting points for the best fit
    defaults = {'lM'   : 14.37+np.log10(h), #Result of SV relation
                'conc'    : 5.0, #Arbitrary
                'tau' : 0.153, #Y1
                'fmis' : 0.32, #Y1
                'Am'    : 1.02, #Y1 approx.
                'B0'   : 0.07, #Y1
                'Rs'   : 2.49,  #Y1; Mpc physical
                'sig_b': 0.005} #Y1 boost scatter
    return defaults
