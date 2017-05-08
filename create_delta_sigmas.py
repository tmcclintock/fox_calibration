"""
Use the Build-Delta-Sigma repository to turn my HMCF curves into delta sigma profiles.
"""
import os, sys
sys.path.insert("../Build-Delta-Sigma/src/wrapper/"
import Build_Delta_Sigma
import matplotlib.pyplot as plt
import numpy as np

inpath = ""

#The radial points of all data
nbins = 50
bins = np.logspace(np.log10(0.01), np.log10(150.0), nbins+1)
R = (bins[:-1]+bins[1:])/2.

#Snapshot characteristics
inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
c = np.linspace(1.0, 0.0, len(linds))
cmaps = ['Reds','Oranges','Greens', 'Blues']

#This is the fox sim cosmology
cosmo = {"h":0.670435,"om":0.31834,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]
