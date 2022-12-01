import sys
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm as norm 
from scipy.stats import genextreme as gev
from scipy.stats import gumbel_l as gum
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.utils as ut
import processing.lens as lens

lens2_gmst_full = gmst.get_gmst_annual_lens2_ensmean()
lens2_tmax_full = lens.get_LENS2_jan_tmax_QNW()

lens2_gmst = lens2_gmst_full.sel(time=slice('1850', '2021'))
lens2_tmax = lens2_tmax_full.sel(time=slice('1850', '2021'))

lens2_gmst_arr = np.tile(lens2_gmst.values, lens2_tmax.shape[0])
lens2_tmax_arr = np.ravel(lens2_tmax.values)

data = np.ravel(lens2_tmax.sel(time=slice('1850', '1880')).values)

cbins = np.arange(22, 33.5, 0.2)
xx = np.linspace(22, 33.5, 100)
hist_2017, bins = np.histogram(data, bins=cbins, density=True)
width = 1.0 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

fig = plt.figure(figsize=(13,8))
plt.bar(center, hist_2017, align='center', width=width, edgecolor='k', facecolor='b', color='blue', alpha = 0.25, label='2017 LENS2')
plt.plot(xx, norm.pdf(xx, *norm.fit(data)), color='k', label='norm')
plt.plot(xx, gev.pdf(xx, *gev.fit(data)), color='b', label='gev')
plt.plot(xx, gum.pdf(xx, *gum.fit(data)), color='r', label='gev')
plt.show()


xx = np.linspace(27, 34, 1000)
normfit = norm.fit(data)
y_norm = 1/norm.sf(xx, *normfit)
gevfit = gev.fit(data)
y_gev = 1/gev.sf(xx, *gevfit)
gumrfit = gum.fit(data)
y_gum = 1/gum.sf(xx, *gumrfit)

u, tau = ut.get_tau(data)

# CI
nboot = 1000
bspreds = np.zeros((nboot, y_norm.size))

for i in range(nboot):
    z = bootstrap(data, n_samples=10)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/y_norm, *normfit_)

yinf, ysup = np.quantile(bspreds, [0.025, 0.975], axis = 0)

fig = plt.figure(figsize=(12,6))
plt.gca().fill_between(y_norm, ysup, yinf, alpha=.25, label='95% confidence interval')
plt.plot(y_norm, xx, color='k', lw=0.8, alpha = 1, label = 'Parametric return period (norm)')
plt.plot(y_gev, xx, color='g', lw=0.8, alpha = 1, label = 'parametric return period (GEV)')
plt.plot(y_gum, xx, color='r', lw=0.8, alpha = 1, label = 'Parametric return period (Gumbel)')
plt.scatter(tau, u, marker='o', facecolor='lightskyblue', edgecolor='blue', color='blue', alpha = 1, label = 'Non parametric return period')
# plot the max value line
plt.axhline(np.max(data), lw=1, color='grey', ls='dotted', label='max value')
plt.grid(lw=0.2, ls='--', color='grey')
plt.legend()
plt.gca().set_xscale('log')
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
plt.xlim([0.9,11000])
plt.ylim([27, 34])
plt.xlabel('Return period (years)')
plt.ylabel('January Tmax (ºC)')
plt.plot()



