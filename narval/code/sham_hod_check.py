# Compare model HOD to actual for Uchuu SHAM mock
# v0.1.0, 2022-03-30 - Code started with snippets from PHYS 788 A3

# Note: Requires first calculating the halo mass function using
# uchuu_mocks/code/calc_hmf.py

# Imports
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf
import sys

# from gp_Mmin_prediction import *


# Define functions for the central and satellite occupation models
def cen_occ_func(logM, fmax, logM_min, sigma_logM):
    return fmax / 2. * (1. + erf((logM - logM_min) / sigma_logM))


def sat_occ_func(logM, logM_sat, alpha, logM_cut, fmax, logM_min, sigma_logM):
    return ((10.**logM / 10**logM_sat)**alpha *
            np.exp(-10.**logM_cut / 10.**logM) *
            cen_occ_func(logM, fmax, logM_min, sigma_logM) / fmax)


start_time = dt.datetime.now()
print("Starting at", start_time)

print("Loading galaxy mock, elapsed time", dt.datetime.now() - start_time)

# Load galaxy mock catalogue
galaxy_mock = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/mocks/uchuu/"
                         "vpeak_scat0.1_z0.70_satmass.dat")

print("Loading Uchuu HMF, elapsed time", dt.datetime.now() - start_time)

# Load galaxy mock catalogue
uchuu_hmf = np.loadtxt("/home/mj3chapm/RSD/uchuu_mocks/data/uchuu_hmf.dat")

# Arbitrary mass bins
mass_bins = np.arange(11., 15.6, 0.1)
mass_bin_centres = np.arange(11.05, 15.55, 0.1)

# Bin by mass
galaxy_mock_cen = galaxy_mock[galaxy_mock[:, 3] == -1, :]
galaxy_mock_sat = galaxy_mock[galaxy_mock[:, 3] != -1, :]
print("N_cen=", galaxy_mock_cen.shape[0])
print("N_sat=", galaxy_mock_sat.shape[0])
print("f_sat=", galaxy_mock_sat.shape[0] / (galaxy_mock_cen.shape[0] +
                                            galaxy_mock_sat.shape[0]))

gal_cen_occ = np.histogram(np.log10(galaxy_mock_cen[:, 1]), bins=mass_bins)[0]
gal_sat_occ = np.histogram(np.log10(galaxy_mock_sat[:, -1]),
                           bins=mass_bins)[0]

print("Loading models, elapsed time", dt.datetime.now() - start_time)

pars = np.loadtxt("/home/mj3chapm/RSD/fit/output/data_products/"
                  "sham_hod_subsample.dat")

print("Starting model calculations, elapsed time",
      dt.datetime.now() - start_time)

# For each model plot the central, satellite and total occupation
fig, axes = plt.subplots(3, 1, figsize=(3.5, 8), sharex=True,
                         gridspec_kw={'hspace': 0}, dpi=300)
for i in range(pars.shape[0]):
    # logM_min = Prediction(pars[i, :])
    logM_min = 13.7
    print("For sample {}, logM_min={}".format(i, logM_min))

    model_cen_occ = cen_occ_func(mass_bin_centres, pars[i, 15],
                                 logM_min, pars[i, 10])
    model_sat_occ = sat_occ_func(mass_bin_centres, pars[i, 7],
                                 pars[i, 8], pars[i, 9],
                                 pars[i, 15], logM_min,
                                 pars[i, 10])

    if i == 0:
        axes[0].plot(mass_bin_centres, model_cen_occ, color="C0",
                     linestyle="-", alpha=0.2, label="Posterior Models")

    else:
        axes[0].plot(mass_bin_centres, model_cen_occ, color="C0",
                     linestyle="-", alpha=0.2)

    axes[1].plot(mass_bin_centres, model_sat_occ, color="C0", linestyle="-",
                 alpha=0.2)
    axes[2].plot(mass_bin_centres, model_cen_occ + model_sat_occ, color="C0",
                 linestyle="-", alpha=0.2)

    if (i + 1) % 10 == 0:
        print("Finished {} samples, elapsed time".format(i+1),
              dt.datetime.now() - start_time)

axes[0].plot(mass_bin_centres, gal_cen_occ / uchuu_hmf[:, 1], color="k",
             linestyle="-", label="Uchuu SHAM Mock")
axes[0].legend()
axes[0].set_ylabel(r"$\langle N_{cen} \rangle$")
axes[0].set_yscale('log')
axes[0].set_ylim(10**-5, 10**1.5)

axes[1].plot(mass_bin_centres, gal_sat_occ / uchuu_hmf[:, 1], color="k",
             linestyle="-")
axes[1].set_ylabel(r"$\langle N_{sat} \rangle$")
axes[1].set_yscale('log')
axes[1].set_ylim(10**-5, 10**1.5)

axes[2].plot(mass_bin_centres, (gal_sat_occ + gal_cen_occ) / uchuu_hmf[:, 1],
             color="k", linestyle="-")
axes[2].set_xlabel(r"$M\,[\log(M_\odot / h)]$")
axes[2].set_ylabel(r"$\langle N_{tot} \rangle$")
axes[2].set_yscale('log')
axes[2].set_xlim(11, 15.5)
axes[2].set_ylim(10**-5, 10**1.5)

plt.tight_layout()
plt.savefig("/home/mj3chapm/RSD/fit/output/plots/"
            "hod_check_em_Mmin.jpg")

print("Script complete! Elapsed time".format(i+1),
      dt.datetime.now() - start_time)
