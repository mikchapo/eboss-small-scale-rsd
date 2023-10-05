# Analyze chains from Abacus emulator fits
# v1.0.1, 2022-06-29 - Updated 'aemulus' to emulator

from collections import OrderedDict
from getdist.chains import ParamError
from getdist.mcsamples import loadMCSamples
import getdist.plots as gdplt
import numpy as np
import os
import seaborn as sns
import sys
import yaml

print("Starting...")

def omegam(z, omegam0=0.315):
    return omegam0 * (z + 1.)**3. / (omegam0 * (z + 1.)**3. + 1. - omegam0)


def fsigma8_approximate(z, sigma8=0.811, omegam0=0.315):
    l = 0.83 + 0.19 / (omegam0**0.8)
    beta = 0.98 + 0.01 / (omegam0**1.2)
    gamma = 0.64 + 0.09 / (omegam0**0.4)
    print("lambda:", l, "beta:", beta, "gamma:", gamma)
    return l * sigma8 * omegam(z, omegam0=omegam0)**(gamma) / (1. + z)**beta


fsig8_ref = fsigma8_approximate(0.7368961822796876)
# paths = ["/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/aemulus_fmax/eboss_ao_large_only_v3/eboss_ao_large_only_v3",
#          "/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/aemulus_fmax/eboss_ao_rm_small_v3/eboss_ao_rm_small_v3",
#          "/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/aemulus_fmax/eboss_ao_rm_large_v3/eboss_ao_rm_large_v3",
#          "/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/aemulus_fmax/eboss_ao_tp_v2/eboss_ao_tp_v2"]
# legend_labels = ["$7-60\,h^{-1}$ Mpc", "$0.8-60\,h^{-1}$ Mpc",
#                  "$0.1-7\,h^{-1}$ Mpc", "$0.1-60\,h^{-1}$ Mpc"]
# colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
paths = ["/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/aemulus_fmax/eboss_ao_rm_large_v3/eboss_ao_rm_large_v3",
         "/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/aemulus_fmax/eboss_ao_rm_small_v3/eboss_ao_rm_small_v3",
         "/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/aemulus_fmax/eboss_ao_large_only_v3/eboss_ao_large_only_v3"]
legend_labels = [r"$0.1-7\,h^{-1}$ Mpc", r"$0.8-60\,h^{-1}$ Mpc",
                 r"$7-60\,h^{-1}$ Mpc"]
colors = ["tab:blue", "tab:orange", "tab:green"]

print("fsig8_ref:", fsig8_ref)

params = ["fsigma8_comp"]
samples = []
for a, path in enumerate(paths):
    folder, name = os.path.split(os.path.abspath(path))
    sample = loadMCSamples(path, settings={"ignore_rows": 15000},
                           no_cache=True)
    samples.append(sample)

fig_width = 5.6

g = gdplt.get_single_plotter(width_inch=fig_width)
g.plot_1d(samples, "fsigma8_comp", normalized=True, colors=colors,
          label=r"$f\sigma_8$")
g.add_x_marker(fsig8_ref)
g.add_legend(legend_labels, legend_loc='upper left')
g.export("scale_dependence.jpg")



# Change Log
# v1.0.0, 2022-01-14 - Copied from RSD/fit/code/analyze_afmax.py
