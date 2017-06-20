"""
This is a file that combines all the small scripts I used earlier.
This is being written so that I can use general paths to the catalogs of 
interest isntead of hard-coding in paths to certain catalogs.

Note: the 'richness calibration' is in old_files/ and the 'mass calibration'
is in old_files/SV_version.
"""
import os, sys
import numpy as np
sys.path.insert(0, "./src/")
from CF_functions import *

inds = [6,7,8,9]
zs = [1.0, 0.5, 0.25, 0.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]
linds = range(0,7)
halobase="/calvin1/tmcclintock/fox_data/richness_halos/rich_snapdir_ps%d_%03d/"
mpath = halobase+"/mass_halos_ps%d_m%d_%03d.txt" #Path to mass split halos
lpath = halobase+"/richness_halos_ps%d_l%d_%03d.txt" #Path to lam split halos
lam_edges = [5, 10, 14, 20, 30, 45, 60, np.inf]
lM_edges = [13.1, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.5, 15.0]#, 16.0]

dmpath = "/calvin1/tmcclintock/down_sampled_snapshots/snapdir_%03d/snapshot_%03d_z%s_down10000"

if __name__ == "__main__":
    pscatter = 25
    calc_hhcf(lpath%(pscatter,6,pscatter,0,6))
    print "HHCF done"
    RR, DdRd, Dd, Rh = calc_DdRd_and_RR(dmpath%(6,6,"1.0"))
    print "DdRd and RR done"
    calc_hmcf(lpath%(pscatter,6,pscatter,0,6), RR, DdRd, Dd, Rh, "test.txt")
    print "All done"

