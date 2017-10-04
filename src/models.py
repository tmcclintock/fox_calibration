"""
This file contains interfaces for the boost factor model and the
DeltaSigma model.
"""
import numpy as np
import os, sys
import clusterwl
import helper_functions as HF
cosmo = HF.get_cosmo()
h = cosmo['h']
om = cosmo['om']
defaults = HF.get_model_defaults(cosmo['h'])

#Swap between whatever model type we are working with and return
#the parameters, including their default values.
def model_swap(params, name):
    c, tau, fmis, Am, B0, Rs, sigb = [defaults['conc'], defaults['tau'], defaults['fmis'], defaults['Am'], defaults['B0'], defaults['Rs'], defaults['sig_b']]
    if name is "full_calibration":
        lM, c, tau, fmis, Am = params
    if name is "Mc":
        lM, c = params
    return [lM, c, tau, fmis, Am, B0, Rs, sigb]

def get_delta_sigma(params, z, Rlam, cosmo, k, Plin, Pnl, Rmodel, xi_mm, Redges, model_name):
    lM, c, tau, fmis, Am, B0, Rs, sigb = model_swap(params, model_name)
    M = 10**lM #Msun/h
    xi_nfw   = clusterwl.xi.xi_nfw_at_R(Rmodel, M, c, om)
    bias = clusterwl.bias.bias_at_M(M, k, Plin, om)
    xi_2halo = clusterwl.xi.xi_2halo(bias, xi_mm)
    xi_hm    = clusterwl.xi.xi_hm(xi_nfw, xi_2halo)
    Rp = np.logspace(-2, 2.4, 1000, base=10) #Mpc/h
    Sigma  = clusterwl.deltasigma.Sigma_at_R(Rp, Rmodel, xi_hm, M, c, om)
    DeltaSigma = clusterwl.deltasigma.DeltaSigma_at_R(Rp, Rp, Sigma, M, c, om)
    Rmis = tau*Rlam #Mpc/h
    Sigma_mis  = clusterwl.miscentering.Sigma_mis_at_R(Rp, Rp, Sigma, M, c, om, Rmis, kernel="exponential")
    DeltaSigma_mis = clusterwl.miscentering.DeltaSigma_mis_at_R(Rp, Rp, Sigma_mis)
    full_profile = (1-fmis)*DeltaSigma + fmis*DeltaSigma_mis
    full_profile *= Am #multiplicative bias

    ave_profile = np.zeros((len(Redges)-1))
    #Note: Redges already in Mpc/h comoving here
    clusterwl.averaging.average_profile_in_bins(Redges, Rp, full_profile, ave_profile)
    return Rp, DeltaSigma, ave_profile
