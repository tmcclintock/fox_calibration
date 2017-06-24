"""
Here I actually do the calibration analysis, calling the optimizer/MCMC for each redshift, richness, and scatter value.
"""
import os, sys
sys.path.insert(0,"./src/")
from likelihoods import *
import numpy as np
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op

if __name__ == "__main__":
    print "not implemented"
