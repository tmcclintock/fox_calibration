import numpy as np
import os, sys
from models import *
import helper_functions as HF
cosmo = HF.get_cosmo()
h = cosmo['h']
om = cosmo['om']

"""
Log prior
Note: depending on what model we have, many of these factors will just be turned off.
"""
def lnprior(params, name):
    lM, c, tau, fmis, Am, B0, Rs, sigb = model_swap(params, name)
    if lM < 11.0 or lM > 18.0 or c <= 0.0 or c > 20.0 or Am <= 0.0 or tau <= 0.0 or Rs <=0.0 or B0 < 0.0 or sigb < 0.0 or fmis < 0.0 or fmis > 1.0: return -np.inf
    #Note: Rlam is Mpc/h here
    LPfmis = (0.32 - fmis)**2/0.05**2 #Y1
    LPtau  = (0.153 - tau)**2/0.03**2 #Y1
    LPA    = (1.02 - Am)**2/0.038**2 #SV
    return -0.5*(LPfmis + LPtau + LPA)

"""
Log posterior of the DeltaSigma model
"""
def lnlike(params, args):
    z, lam, Rlam, Rdata, ds, icov, cov, cosmo, k, Plin, Pnl, Rmodel, xi_mm, Redges, indices, model_name = args
    lM, c, tau, fmis, Am, B0, Rs, sigb = model_swap(params, model_name)

    LLDS = 0
    Rp, full_DeltaSigma, ave_DeltaSigma = get_delta_sigma(params, z, Rlam, cosmo, k, Plin, Pnl, Rmodel, xi_mm, Redges, model_name)
    ds_model = ave_DeltaSigma[indices]
    #print indices
    #sys.exit()
    ds_model *= h*(1+z)**2 #make the model physical
    X = ds - ds_model
    LLDS = -0.5*np.dot(X, np.dot(icov, X))
    return LLDS #+ LLboost #no boost factors in fake data


"""
Log posterior probability
"""
def lnprob(params, args):
    z, lam, Rlam, Rdata, ds, icov, cov, cosmo, k, Plin, Pnl, Rmodel, xi_mm, Redges, indices, model_name = args
    lpr = lnprior(params, model_name)
    if not np.isfinite(lpr): return -np.inf
    return lpr + lnlike(params, args)
