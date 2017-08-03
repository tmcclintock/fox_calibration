import numpy as np
import os, sys, emcee
import scipy.optimize as op
sys.path.insert(0, "../Delta-Sigma/src/wrapper/")
import py_Delta_Sigma as pyDS
from colossus.halo import concentration as conc
from colossus.cosmology import cosmology as col_cosmology

#This is the fox sim cosmology
h = 0.670435
cosmo = {"h":h,"om":0.31834,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]
colcos = {"H0":cosmo['h']*100.,"Om0":cosmo['om'], 
          'Ob0': 0.049017, 'sigma8': 0.83495, 'ns': 0.96191, 'flat':True}
col_cosmology.addCosmology('fiducial_cosmology', colcos)
col_cosmology.setCosmology('fiducial_cosmology')
input_params = {"NR":300,"Rmin":0.01,
                "Rmax":200.0,"Nbins":15,"delta":200,
                "Rmis":0.0, "fmis":0.0,
                "miscentering":0,"averaging":1}

#Write the likelihoods
def lnprior(params):
    lM = params
    #Mass to low or too high
    if lM < 10 or lM > 18: return -np.inf
    return 0.0

def lnlike(params, R, ds, icov, flags, z, extras):
    lM = params
    klin, knl, Plin, Pnl, cosmo, inparams = extras
    inparams['Mass'] = 10**lM #Mpc/h
    inparams["concentration"] = conc.concentration(10**lM, '200m', z, model='diemer15')
    result = pyDS.calc_Delta_Sigma(klin, Plin, knl, Pnl, cosmo, inparams)
    model = result['ave_delta_sigma']*h*(1.+z)**2 #Msun/pc^2 physical
    model = model[flags]
    X = ds - model
    return -0.5*np.dot(X, np.dot(icov, X))

def lnprob(params, R, ds, icov, flags, z, extras):
    lM = params
    lpr = lnprior(params)
    if not np.isfinite(lpr): return -np.inf
    return lpr + lnlike(params, R, ds, icov, flags, z, extras)

