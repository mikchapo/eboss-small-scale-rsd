"""
Analyze chains from SHAM tests of varying number density and redshift.

v0.1.0, 2023-04-26 - Copied from P2/fit/code/analyze_chains.py
"""

from getdist.mcsamples import loadMCSamples
import getdist.plots as gdplt
import numpy as np
import os


path_root = ("/home/mj3chapm/projects/rrg-wperciva/mj3chapm/RSD/chains/"
             "joint_planck_aemulus")
file_mocks = ["z057_n1", "z07_n4", "z07_n1"]
mocks = ["z=0.57, n=1e-4", "z=0.7, n=4e-4", "z=0.7, n=1e-4"]
colors = ["tab:blue", "tab:orange", "tab:green"]

samples = []
for i in range(len(file_mocks)):
    run_name = "uchuu_{}_scale_s8_v5".format(file_mocks[i])
    path = "{}/{}/jpa_{}".format(path_root, run_name, run_name)
    folder, name = os.path.split(os.path.abspath(path))
    sample = loadMCSamples(path, settings={"ignore_rows": 0.3},
                           no_cache=True)
    p = sample.getParams()
    sample.addDerived(p.fsigma8_comp, name='fsigma8_new', label=r'f\sigma_{8}}')
    samples.append(sample)

scaling = np.sqrt(3)
width = 6.375

gdplot = gdplt.get_single_plotter(width_inch=width*scaling)
gdplot.settings.axes_fontsize = 10*scaling
gdplot.settings.axes_labelsize = 10*scaling
gdplot.settings.legend_fontsize = 10*scaling
gdplot.plot_1d(samples, "fsigma8_new", colors=colors)
gdplot.add_x_marker(0.4726383851649886)
gdplot.add_legend(mocks)
gdplot.export("../output/plots/sham_n-z_test_scaled.png")

# Change Log
