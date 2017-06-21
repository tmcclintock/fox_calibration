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

#The radial locations of the actual data.
nbins = 15
bins = np.logspace(np.log10(0.0323), np.log10(30.0), nbins+1)
R_data = (bins[:-1]+bins[1:])/2. #Locations of the data bins; Mpc physical

def calc_DS(R, xi_hm, Mass, redshift, savepath=None, cosmo=None):
    if not cosmo:
        #This is the fox sim cosmology
        cosmo = {"h":0.670435,"om":0.31834,"ok":0.0}
        cosmo["ode"]=1.0-cosmo["om"]
        #Here is the same cosmology but for cosmocalc
        colcos = {"H0":cosmo['h']*100.,"Om0":cosmo['om'], 
                  'Ob0': 0.049017, 'sigma8': 0.83495, 'ns': 0.96191, 'flat':True}
        col_cosmology.addCosmology('fiducial_cosmology', colcos)
        col_cosmology.setCosmology('fiducial_cosmology')

    params = {"Mass": Mass, "delta":200, 
              "timing":1, "miscentering":0}
    params["concentration"] = conc.concentration(Mass, '200m', redshift, 
                                                 model='diemer15')
    results = BDS.build_Delta_Sigma(R, xi_hm, cosmo, params)
    DS = results['delta_sigma']
    if savepath:
        out = np.array([R,DS]).T
        np.savetxt(savepath, out)
    return [R, DS]

def create_data_vector(Rinit, DSinit, C, z, save, Csave, cosmo=None):
    if 1 == 1 : #developing still
        return
    if not cosmo:
        #This is the fox sim cosmology
        cosmo = {"h":0.670435,"om":0.31834,"ok":0.0}
        cosmo["ode"]=1.0-cosmo["om"]
    h = cosmo['h']
    #change units to Msun,Mpc physical
    R = Rinit/(h*(1+z))
    DS = DSinit*h*(1+z)**2
    #Spline to get the curve at R_data
    sp = IUS(R,DS)
    dsout = spl(R_data)
    np.savetxt(save, np.array(R,dsout))
    np.savetxt(Csave, C)
    return
