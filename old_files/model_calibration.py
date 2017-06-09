"""
Model the calibration.
"""
import numpy as np
import emcee, sys
import scipy.optimize as op

cut = 3 #don't use any measurements below this
M_true = np.genfromtxt("txt_files/true_masses.txt")[:,cut:]
lam = np.genfromtxt("txt_files/true_richnesses.txt")[:,cut:]
lM = np.genfromtxt("txt_files/mcmc_masses.txt")[:,cut:]
lMe = np.genfromtxt("txt_files/mcmc_errs.txt")[:,cut:]
C = M_true/10**lM
Ce = np.log(10) * lMe * C
zs = [1.0, 0.5, 0.25, 0.0]
z = np.copy(zs)
for i in range(3):
    z = np.vstack((z, zs))
z = np.array(z).T
print z.shape
print lam.shape
print C.shape

z0p1 = 1.5 #zpivot = 0.5
h = 0.670435 #hubble constant/100
Mpivot = 10**(13.8 + np.log10(h))
lampivot = 30.0
print Mpivot

def get_C_model(params, lam, z):
    C0, a, b = params
    return C0*((1.+z)/z0p1)**a*(lam/lampivot)**b

def lnprior(params):
    C0, a, b = params
    if C0 < 0.0: return -np.inf
    return 0.0

def lnlike(params, C, Cerr, lam, z):
    Cmodel = get_C_model(params, lam, z)
    return -0.5*np.sum((C-Cmodel)**2/Cerr**2)

def lnprob(params, C, Cerr, lam, z):
    lnp = lnprior(params)
    if not np.isfinite(lnp): return -np.inf
    return lnp + lnlike(params, C, Cerr, lam, z)

lnprobargs = (C, Ce, lam, z)
nll = lambda *args: -lnprob(*args)
result = op.minimize(nll, [1.05, 0.2, 0.2], args=lnprobargs)
print result

import emcee
nwalkers = 8
ndim = 3
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=lnprobargs)
nsteps = 10000
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler.run_mcmc(pos, nsteps)
chain = sampler.flatchain
np.savetxt("txt_files/calchain.txt", chain)
means = np.mean(chain, 0)
stds = np.std(chain, 0)
print means, stds

