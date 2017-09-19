import os, sys
sys.path.insert(0,"./src/")
from get_data import *
sys.path.insert(0, "../Delta-Sigma/src/wrapper/")
import py_Delta_Sigma as pyDS
from colossus.halo import concentration as conc
from colossus.cosmology import cosmology as col_cosmology
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", size=20)

#Fox cosmology
h = 0.670435
cosmo, input_params = get_cosmo_and_params()
colcos = get_colcos()
col_cosmology.addCosmology('fiducial_cosmology', colcos)
col_cosmology.setCosmology('fiducial_cosmology')

def plot_data(ps, zi, lj, use_y1=True):
    z = zs[zi]
    klin, Plin, knl, Pmm = get_P(zi)
    input_params["R_bin_min"] = 0.0323*(h*(1+z)) #Mpc/h comoving
    input_params["R_bin_max"] = 30.0*(h*(1+z)) #Mpc/h comoving
    extras = [klin, knl, Plin, Pmm, cosmo, input_params]
    R, DS, cov, cut = get_data_and_cov(ps, zi, lj, use_y1, nocut=True)
    DSe = np.sqrt(np.diag(cov))
    bad = R < 0.2 #Mpc
    gud = R > 0.2 #Mpc
    if use_y1: c = 'b'
    else: c = 'g'
    plt.errorbar(R[bad], DS[bad], DSe[bad], ls='', marker='.', mfc='w', c=c)
    plt.errorbar(R[gud], DS[gud], DSe[gud], ls='', marker='.', c=c)
    return R, DS, DSe

def plot_model(M, zi):
    z = zs[zi]
    inparams = input_params.copy()
    inparams['Mass'] = M #Msun/h
    inparams["concentration"] = conc.concentration(M, '200m', z, model='diemer15')
    klin, Plin, knl, Pmm = get_P(zi)
    input_params["R_bin_min"] = 0.0323*(h*(1+z)) #Mpc/h comoving
    input_params["R_bin_max"] = 30.0*(h*(1+z)) #Mpc/h comoving
    extras = [klin, knl, Plin, Pmm, cosmo, input_params]
    result = pyDS.calc_Delta_Sigma(klin, Plin, knl, Pmm, cosmo, inparams)
    R = result["Rbins"]/(h*(1+z)) #Mpc physical
    DS = result['ave_delta_sigma']*h*(1+z)**2 #Msun/pc^2 physical
    return R, DS

def plot_bf(ps, zi, lj, use_y1=True):
    bfM = get_bf_results(ps, use_y1)
    R, DS = plot_model(bfM[zi, lj], zi)
    plt.loglog(R, DS, c='r', label="Best", ls='--')
    return R, DS

def plot_truth(ps, zi, lj):
    trueM = get_true_masses(ps)
    plot_model(trueM[zi, lj], zi)
    R, DS = plot_model(trueM[zi, lj], zi)
    plt.loglog(R, DS, c='blue', label="Truth")
    return R, DS

def print_info(ps, zi, lj, use_y1=True):
    trueM = get_true_masses(ps)
    bfM = get_bf_results(ps, use_y1)
    bfcal = trueM/bfM
    print "Best fit mass = ", np.log10(bfM[zi, lj])
    print "Truth mass = ", np.log10(trueM[zi, lj])
    print "BF calibration = ", bfcal[zi, lj]

def get_title(ps, zi, lj, use_y1=True):
    title = r"${\rm scatter=\ }%d\ z=%.2f\ l%d\ C_{\rm %s}$"
    if use_y1: s = "Y1"
    else: s = "SV"
    return title%(ps, zs[zi], lj, s)

def plot_comparison(R, DS, DSe, Rm, DSm, c='k'):
    from scipy.interpolate import InterpolatedUnivariateSpline as IUS
    model = IUS(Rm, DSm)
    pdiff = (DS-model(R))/model(R)
    epd   = DSe/model(R)
    bad = R < 0.2 #Mpc
    gud = R > 0.2 #Mpc
    plt.errorbar(R[bad], pdiff[bad], epd[bad], ls='', marker='.', mfc='w', c=c)
    plt.errorbar(R[gud], pdiff[gud], epd[gud], ls='', marker='.', c=c)
    plt.axhline(0, c='k')
    return

if __name__ == "__main__":
    use_y1 = True
    fig, axarr = plt.subplots(2, sharex=True)
    plt.sca(axarr[0])
    R, DS, DSe = plot_data(0, 0, 4, False)
    R, DS, DSe = plot_data(0, 0, 4, use_y1)
    Rt, DSt = plot_truth(0, 0, 4)
    Rbf, DSbf = plot_bf(0, 0, 4, use_y1)
    plt.legend(loc=0, fontsize=12)
    plt.sca(axarr[1])
    plot_comparison(R, DS, DSe, Rt, DSt, c="blue")
    plot_comparison(R, DS, DSe, Rbf, DSbf, c="red")
    print_info(0, 0, 4, use_y1)
    title = get_title(0, 0, 4, use_y1)
    #plt.title(title)
    plt.xscale('log')
    plt.yscale('log')
    axarr[1].set_yscale('linear')
    axarr[1].set_xlabel(r"$R\ [{\rm Mpc}]$")
    axarr[1].set_ylabel("Frac. Diff.")
    axarr[0].set_ylabel(r"$\Delta\Sigma\ [{\rm M_\odot/ pc^2}]$")
    plt.subplots_adjust(bottom=0.15, left=0.15, hspace=0.0)
    plt.savefig("y1sv_comparison.png"
    plt.show()
