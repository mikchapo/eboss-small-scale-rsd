"""
Plot a comparison of DES Y1 contours and eBOSS LRG S8 constraints.

v0.3.0, 2023-07-17 - Code copied from Graham and updated for thesis
"""


from getdist.mcsamples import loadMCSamples
import copy
import getdist.plots as gdplt
import numpy as np
import os
import sys

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

# output = sys.argv[1]
# ignore_rows = float(sys.argv[2])
# paths = sys.argv[3:]

output = "/home/mj3chapm/RSD/fit/output/plots/thesis/afmax_eboss_des_comp"
ignore_rows = 20000.
paths = ["/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/des/des_default",
         "/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/aemulus_fmax/eboss_ao_tp_v2/eboss_ao_tp_v2",
         "/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/aemulus_fmax/eboss_ao_large_only_v3/eboss_ao_large_only_v3",
         "/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/planck_test_2/planck_test"]

folders = []
names = []
full_fit_samples = []
quasi_lin_samples = []
for path in paths:
    folder, name = os.path.split(os.path.abspath(path))
    folders.append(folder)
    names.append(name)
    if name == "planck_test":
        sample = loadMCSamples(path)       

    elif name == "des_default":
        sample = loadMCSamples(path, settings={"ignore_rows": 30000}, no_cache=True)

    else:
        sample = loadMCSamples(path, settings={"ignore_rows": ignore_rows}, no_cache=True)



    if name == "des_default" or name == "planck_test":
        p = sample.getParams()
        S8 = p.s8omegamp5 / (0.3**0.5)
        sample.addDerived(S8, name='S8', label=r'S_8')
        full_fit_samples.append(sample)
        quasi_lin_samples.append(sample)

    else:
        p = sample.getParams()
        S8 = p.sigma8 * np.power(p.omegam/0.3, 0.5)
        scaled_sample = copy.deepcopy(sample)
        S8_scaled = p.f * S8
        sample.addDerived(S8, name='S8', label=r'S_8')
        scaled_sample.addDerived(S8_scaled, name='S8', label=r'S_8')
        if name == "eboss_ao_tp_v2":
            full_fit_samples.append(scaled_sample)
            full_fit_samples.append(sample)
        else:
            quasi_lin_samples.append(scaled_sample)
            quasi_lin_samples.append(sample)

    mean = sample.mean("S8")
    var = sample.var("S8")
    print(name)
    print("S8", mean, var**0.5)
    print("\n")


gdplot_settings = gdplt.GetDistPlotSettings(fig_width_inch=3.175)
# gdplot = gdplt.get_single_plotter(settings=gdplot_settings)
# gdplot.plot_1d(samples, 'S8')
# gdplot.add_legend(legend_labels=["DES Y1", "Small-scale eBOSS LRG", "Planck18"])
# gdplot.export(output + "_1d.jpg")

gdplot = gdplt.get_subplot_plotter(settings=gdplot_settings)
gdplot.triangle_plot(full_fit_samples, ["omegam", "S8"], filled=True, legend_labels=["DES Y1", r"$\gamma_f S_8$ eBOSS Full", r"eBOSS Full", "Planck18"],
                     contour_colors=["tab:orange", "tab:green", "tab:blue", "tab:red"])
gdplot.export(output + "_full_fit.jpg")

gdplot = gdplt.get_subplot_plotter(settings=gdplot_settings)
gdplot.triangle_plot(quasi_lin_samples, ["omegam", "S8"], filled=True, legend_labels=["DES Y1", r"$\gamma_f S_8$ eBOSS Q.-L.", r"eBOSS Q.-L.", "Planck18"],
                     contour_colors=["tab:orange", "tab:green", "tab:blue", "tab:red"])
gdplot.export(output + "_quasi_lin.jpg")

# gdplot = gdplt.get_subplot_plotter(settings=gdplot_settings)
# gdplot.triangle_plot(samples, ["omegam", "S8_f"], filled=True, legend_labels=["DES Y1", "Small-scale eBOSS LRG", "Planck18"], contour_colors=["tab:orange", "tab:blue", "tab:green"])
# gdplot.export(output + "_S8_f_2d.jpg")
