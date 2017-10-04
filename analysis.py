"""
Here I actually do the calibration analysis, calling the optimizer/MCMC for each redshift, richness, and scatter value.
"""
import os, sys
sys.path.insert(0,"./src/")
from likelihoods import *
from helper_functions import *
import numpy as np
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op
import clusterwl

use_y1 = False

#MCMC parameters
nwalkers = 4
ndim = 1
nsteps = 4000
nburn = nsteps/4

#Fox cosmology
cosmo = get_cosmo()
h = cosmo['h']
om = cosmo['om']
zs, zstrings = get_zs()

def do_best_fit(bf_args, ps, zi, lj, true_M):
    z, lam, Rlam, Rdata, ds, icov, cov, cosmo, k, Plin, Pnl, Rmodel, xi_mm, Redges, indices, model_name = bf_args
    true_lM = np.log10(true_M)

    if model_name is 'Mc':guess = [true_lM, 4.1] 
    elif model_name is "full_calibration": guess = [true_lM, 4.1, 0.153, 0.32, 1.02]
    nll = lambda *args: -lnprob(*args)
    result = op.minimize(nll, x0=guess, args=bf_args)
    print result
    np.savetxt("output_files/bf_files/bf_ps%d_z%d_l%d.txt"%(ps, zi, lj), result['x'])
    return [10**result['x'][0], result['x'][1]]

def do_mcmc():
    """
                sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(R, DS, icov, z, extras), threads=4)
                pos = [bf_masses[i, j] + 1e-4*np.random.randn(ndim) for nw in range(nwalkers)]
                print "Running chain for z=%f l%d starting at lM=%.3f"%(z, j, bf_masses[i,j])
                sampler.run_mcmc(pos, nsteps)
                print "\tchain complete for z=%f l%d starting at lM=%.3f"%(z, j, bf_masses[i,j])
                fullchain = sampler.flatchain
                likes = sampler.flatlnprobability
                np.savetxt("output_files/chains/chain_ps%d_z%.2f_l%d.txt"%(ps, z, j), fullchain)
                np.savetxt("output_files/chains/likes_ps%d_z%.2f_l%d.txt"%(ps, z, j), likes)
                chain = fullchain[nburn*nwalkers:]
                mcmc_masses[i,j] = np.mean(chain)
                mcmc_stds[i,j] = np.std(chain)
                cal[i,j] = 10**true_lM[i,j]/10**mcmc_masses[i,j]
                calerr[i,j] = np.log(10) * mcmc_stds[i,j] * cal[i,j]
                print "MCMC done for ps%d z%d, l%d"%(ps, i, j)
                print "mcmc result is ",mcmc_masses[i,j], cal[i,j]
                continue #end j
            continue #end i,ind
        np.savetxt("output_files/mass_fits/mcmc_masses_ps%d.txt"%ps, mcmc_masses)
        np.savetxt("output_files/mass_fits/mcmc_cal_ps%d.txt"%ps, cal)
        np.savetxt("output_files/mass_fits/mcmc_calerr_ps%d.txt"%ps, calerr)
        continue #end ps
        """
    return 0

if __name__ == "__main__":
    save = True
    for ps in [0]:#[15, 25, 35]:
        lams = get_lams(ps)
        Rlams = (lams/100.0)**0.2 #Mpc/h; richness radius
        true_M = get_true_M(ps)
        bf_M   = np.ones_like(true_M)
        bf_c   = np.ones_like(true_M)
        bf_cal = np.ones_like(true_M)
        for i in xrange(0,4):
            z = zs[i]
            k, Plin, Pnl = get_P(z)
            Rmodel = np.logspace(-2, 3, num=1000, base=10) 
            xi_mm = clusterwl.xi.xi_mm_at_R(Rmodel, k, Pnl)
            for j in xrange(0,7):#linds:
                print "Starting analysis for z%d z=%.2f l%d"%(i, z, j)
                lam = lams[i, j]
                Rlam = Rlams[i, j]
                #Xi_mm MUST be evaluated to higher than BAO for correct accuracy
                model_name = "full_calibration"

                Rdata, ds, icov, cov, inds = get_data_and_cov(ps, i, j, use_y1=use_y1, nocut=False)
                Redges = get_Redges(use_y1)*h*(1+z) #Made Mpc/h comoving
                args = [z, lam, Rlam, Rdata, ds, icov, cov, cosmo, k, Plin, Pnl, Rmodel, xi_mm, Redges, inds, model_name]
                
                bf_M[i, j], bf_c[i, j] = do_best_fit(args, ps, i, j, true_M[i,j])
                bf_cal[i,j] = true_M[i, j]/bf_M[i, j]
                print "Best fit done for ps%d z%d z=%.2f l%d"%(ps, i, z, j)
                print "Bf is ",bf_M[i, j], bf_cal[i, j]
                continue #end j
            continue #end i
        if save:
            if use_y1:
                np.savetxt("output_files/mass_fits/bf_y1_masses_ps%d.txt"%ps, bf_M)
                np.savetxt("output_files/mass_fits/bf_y1_conc_ps%d.txt"%ps, bf_c)
                np.savetxt("output_files/mass_fits/bf_y1_cal_ps%d.txt"%ps, bf_cal)
            else:
                np.savetxt("output_files/mass_fits/bf_sv_masses_ps%d.txt"%ps, bf_M)
                np.savetxt("output_files/mass_fits/bf_sv_conc_ps%d.txt"%ps, bf_c)
                np.savetxt("output_files/mass_fits/bf_sv_cal_ps%d.txt"%ps, bf_cal)

    #do_mcmc()
    #see_results() #NEED TO IMPLEMENT
