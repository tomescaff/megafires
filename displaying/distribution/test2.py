import sys
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gennorm as norm
from scipy.stats import genextreme as gev
from scipy.stats import gumbel_r as gum

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.utils as ut
import processing.lens as lens
import processing.math as pmath

ac_year = '2017'
cf_year = '2070'

# raw values

data_raw = xr.open_dataset('../../../../LENS2_tasmax.nc')['tasmax']-273.15
data_raw = data_raw.where(data_raw.time.dt.month ==1, drop=True)

data = np.ravel(data_raw.sel(time=slice('1850', '1880')).values)

cbins = np.arange(22, 36.5, 0.2)
xx = np.linspace(22, 36.5, 100)
hist_2017, bins = np.histogram(data, bins=cbins, density=True)
width = 1.0 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2


fig = plt.figure(figsize=(13,8))
plt.bar(center, hist_2017, align='center', width=width, edgecolor='k', facecolor='b', color='blue', alpha = 0.25, label='2017 LENS2')
plt.plot(xx, norm.pdf(xx, *norm.fit(data)), color='k', label='norm')
plt.plot(xx, gev.pdf(xx, *gev.fit(data)), color='b', label='gev')
plt.plot(xx, gum.pdf(xx, *gum.fit(data)), color='r', label='gev')
plt.show()


xx = np.linspace(27, 36.5, 1000)
normfit = norm.fit(data)
y_norm = 1/norm.sf(xx, *normfit)
gevfit = gev.fit(data)
y_gev = 1/gev.sf(xx, *gevfit)
gumrfit = gum.fit(data)
y_gum = 1/gum.sf(xx, *gumrfit)

u, tau = ut.get_tau(data)


fig = plt.figure(figsize=(12,6))
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
plt.ylim([27, 36.5])
plt.xlabel('Return period (years)')
plt.ylabel('January Tmax (ºC)')
plt.plot()



