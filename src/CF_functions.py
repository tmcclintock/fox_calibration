"""
This contains the functions needed to compute HH and HM correlation
functions. At present it only does HH autocorrelations, not
cross correllations between different halo populations. Also,
no jackknifing is done.
"""
import os, sys
import numpy as np
import pygadgetreader as pgr
from Corrfunc.theory.xi import xi
from Corrfunc.theory.DD import DD
from Corrfunc.utils import convert_3d_counts_to_cf

#Default Corrfunc configuration
config_default = {
    "boxsize": 1050.0,
    "nthreads": 8,
    "nbins": 50}
config_default["bins"] = np.logspace(np.log10(0.01), np.log10(150.0), nbins+1)


def calc_hhcf(halo_path, save_path=None, config=config_default):
    """Calculate the halo-halo correlation function
    Note: assumes the halo catalog has the format X,Y,Z,N,M,Richness.

    Args:
        halo_path (string): Path to the halo catalog
        save_path (string, optional): Path to save the output
        config (dictionary): Configuration of the Corrfunc run. Default is
        config_default
    """
    boxsize  = config["boxsize"]
    nthreads = config["nthreads"]
    bins     = config["bins"]
    X, Y, Z, Np, M, Rich = np.genfromtxt(halo_path, unpack=True)
    result = xi(boxsize, nthreads, bins, X, Y, Z)
    if save_path: np.savetxt(save_path, result)
    return

def calc_DdRd_and_RR(dm_path, N_rand=None, M=1.5, config=config_default):
    """Calculate the DdRd and RR terms for a HMCF. This is because
    these terms take the longest, so it should just be done once
    if the HMCF is done for many mass bins.

    Args:
        dm_path (string): Path to DM catalog
        N_rand (int, optional): number of random points to use. If not specified, then 1.5*N_dm is used.
        M (float): Multiplicative factor to go from N_dm to N_rand. Default is 1.5
        config (dictionary): Configuration of the Corrfunc run. Default is
        config_default

    Returns:
        RR (:obj: Corrfunc.theory.DD): Pair count object from Corrfunc
        DdRd (:obj: Corrfunc.theory.DD): Pair count object from Corrfunc
        Dd (:obj: numpy.array): Numpy array containing the dm catalog
        Rh (:obj: numpy.array): Numpy array containing the halo randoms
    """
    boxsize  = config["boxsize"]
    nthreads = config["nthreads"]
    bins     = config["bins"]
    Xd, Yd, Zd = pgr.readsnap(dmpath%(i, i, zstring), 'pos', 'dm').T
    Xd = Xd.astype('float64')
    Yd = Yd.astype('float64')
    Zd = Zd.astype('float64')
    N_dm = len(Xd)
    if not N_rand: N_rand = int(N_dm*M)
    X1 = np.random.uniform(0, boxsize, N_rand)
    Y1 = np.random.uniform(0, boxsize, N_rand)
    Z1 = np.random.uniform(0, boxsize, N_rand)
    X2 = np.random.uniform(0, boxsize, N_rand)
    Y2 = np.random.uniform(0, boxsize, N_rand)
    Z2 = np.random.uniform(0, boxsize, N_rand)
    RR = DD(0, nthreads, bins, X1, Y1, Z1, X2=X2, Y2=Y2, Z2=Z2)
    DdRd = DD(0, nthreads, bins, Xd, Yd, Zd, X2=X2, Y2=Y2, Z2=Z2)
    return [RR, DdRd, [Xd, Yd, Zd], [X1, Y1, Z1]]

def calc_hmcf(halo_path, RR, DdRd, Dd, Rh, save_path, config=config_default):
    """Calculate the halo-halo correlation function
    Note: assumes the halo catalog has the format X,Y,Z,N,M,Richness.

    Args:
        halo_path (string): Path to the halo catalog
        dm_path (string): Path to the DM catalog
        save_path (string): Path to save the output
        config (dictionary): Configuration of the Corrfunc run. Default is
        config_default
    """
    boxsize  = config["boxsize"]
    nthreads = config["nthreads"]
    bins     = config["bins"]
    Xd, Yd, Zd = Dd
    X1, Y1, Z1 = Rd
    X, Y, Z, Np, M, Rich = np.genfromtxt(halo_path, unpack=True)
    N_h    = len(X)
    N_dm   = len(Xd)
    N_rand = len(X1)
    DhDd = DD(0, nthreads, bins, X1=Xh, Y1=Yh, Z1=Zh, X2=Xd, Y2=Yd, Z2=Zd)
    DhRh = DD(0, nthreads, bins, Xh, Yh, Zh, X2=Xr1, Y2=Yr1, Z2=Zr1)
    hmcf = convert_3d_counts_to_cf(N_h, N_dm, N_rand, N_rand, DhDd, DhRh, DdRd, RR)
    np.savetxt(savepath, hmcf)
    return
