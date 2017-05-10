"""
Use class to get the power spectra for each bin.
We need both P_lin and P_mm (takahashi).
"""
import numpy as np
from classy import Class

#Get the mean redshifts
zs = [1.0, 0.5, 0.25, 0.0]

#Fox cosmology
Ob = 0.049017
Om = 0.31834
Ocdm = Om - Ob
params = {
        'output': 'mPk',
        "h":0.670435,
        "A_s":2.1e-9,
        "n_s":0.96191,
        "omega_b":Ob,
        "omega_cdm":Ocdm,
        'YHe':0.24755048455476272,#By hand, default value
        'P_k_max_h/Mpc':20000.,
        'z_max_pk':1.0,
        'non linear':'halofit'}

cosmo = Class()
cosmo.set(params)
cosmo.compute()


k = np.logspace(-5, 4, base=10, num=1000)
np.savetxt("txt_files/P_files/k.txt", k)
for i in range(len(zs)):
    z = zs[i]

    Pmm  = np.array([cosmo.pk(ki, z) for ki in k])
    Plin = np.array([cosmo.pk_lin(ki, z) for ki in k])
    np.savetxt("txt_files/P_files/Pnl_z%.2f.txt"%(z), Pmm)
    np.savetxt("txt_files/P_files/Plin_z%.2f.txt"%(z), Plin)
    print "Done with z%.2f"%z
