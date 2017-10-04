"""
This contains the functions needed to compute the DeltaSigma curves from
each xi_hm, and then turn those curves into actual data vectors.
"""
import os, sys
import numpy as np
sys.path.insert(0, "../Build-Delta-Sigma/src/wrapper/")
import py_Build_Delta_Sigma as BDS
from colossus.halo import concentration as conc
from colossus.cosmology import cosmology as col_cosmology
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
import clusterwl
import helper_functions as HF
cosmo = HF.get_cosmo()
h = cosmo['h']
om = cosmo['om']

def get_conc(M, z):
    colcos = {"H0":cosmo['h']*100.,"Om0":cosmo['om'], 
              'Ob0': 0.049017, 'sigma8': 0.83495, 'ns': 0.96191, 'flat':True}
    col_cosmology.addCosmology('fox_cosmology', colcos)
    col_cosmology.setCosmology('fox_cosmology')

    om = cosmo['om']
    return conc.concentration(M, '200m', z, model='diemer15')

#The radial locations of the actual data.
def get_binning(y1_binning = True): #Mpc physical
    nbins = 15
    if y1_binning: Redges = np.logspace(np.log10(0.0323), np.log10(30.0), nbins+1)
    else: Redges = np.logspace(np.log10(0.02), np.log10(30.0), nbins+1) #SV binning
    R_mid = (Redges[:-1]+Redges[1:])/2. #Locations of the data bin edges; Mpc physical
    return Redges, R_mid

def calc_DS(R, xi_hm, Mass, z, Rmodel, xi_mm, k, Plin, Pnl, Rlam, dssave=None, avsave=None, y1_binning=True):
    c = get_conc(Mass, z)
    xi_nfw   = clusterwl.xi.xi_nfw_at_R(Rmodel, Mass, c, om)
    bias = clusterwl.bias.bias_at_M(Mass, k, Plin, om)
    xi_2halo = clusterwl.xi.xi_2halo(bias, xi_mm)
    xi_hm_model    = clusterwl.xi.xi_hm(xi_nfw, xi_2halo)

    low  = Rmodel < 0.1 #R[0] #Things are fucked up lower than 0.1 Mpc
    high = Rmodel > R[-1]
    datacut = R > 0.1
    Rall = np.concatenate((Rmodel[low], R[datacut]))
    Rall = np.concatenate((Rall, Rmodel[high]))
    xi_hm_full = np.concatenate((xi_hm_model[low], xi_hm[datacut]))
    xi_hm_full = np.concatenate((xi_hm_full, xi_hm_model[high]))

    Rp = np.logspace(-2, 2.4, 1000, base=10) #Mpc/h
    Sigma  = clusterwl.deltasigma.Sigma_at_R(Rp, Rall, xi_hm_full, Mass, c, om)
    DeltaSigma = clusterwl.deltasigma.DeltaSigma_at_R(Rp, Rp, Sigma, Mass, c, om)
    #Add in miscentering and multiplicative bias
    Am = 1.02 #SV central
    fmis = 0.32 #Y1 central
    tau = 0.153 #Y1 central
    Rmis = tau*Rlam #Mpc/h
    Sigma_mis  = clusterwl.miscentering.Sigma_mis_at_R(Rp, Rp, Sigma, Mass, c, om, Rmis, kernel="exponential")
    DeltaSigma_mis = clusterwl.miscentering.DeltaSigma_mis_at_R(Rp, Rp, Sigma_mis)
    full_profile = (1-fmis)*DeltaSigma + fmis*DeltaSigma_mis
    full_profile *= Am #multiplicative bias
  
    Redges, Rmid = get_binning(y1_binning) #Redges in Mpc physical
    ave_profile = np.zeros((len(Redges)-1))
    #Convert units to Msun, Mpc physical and then calculate the averages
    clusterwl.averaging.average_profile_in_bins(Redges, Rp/(h*(1+z)), full_profile*h*(1+z)**2, ave_profile)

    if dssave:  np.savetxt(dssave, np.array([Rp, full_profile]).T)
    if avsave:  np.savetxt(avsave, np.array([Rmid, ave_profile]).T)
    return [Rp, DeltaSigma, Rmid, ave_profile]
