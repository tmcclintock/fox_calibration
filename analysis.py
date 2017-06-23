"""
Here I actually do the calibration analysis, calling the optimizer/MCMC for each redshift, richness, and scatter value.
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op
