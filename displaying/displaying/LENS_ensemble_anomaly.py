import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys

sys.path.append('../../processing')

import processing.utils as ut

# get LENS time series
basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
lens = da[:, ::12]

# get LENS ensemble anomalies
lens_ensmean = lens.mean('run')
lens_anom = lens - lens_ensmean

# fix LENSA
lens_mean_1920_1940 = lens.sel(time=slice('1920', '1940')).mean('time')
lens_anom_mean_1920_1940 = lens_anom.sel(time=slice('1920', '1940')).mean('time')
lensa = lens_anom - lens_anom_mean_1920_1940 + lens_mean_1920_1940

x = lens_ensmean.time.dt.year.values

fig = plt.figure(figsize=(12,8))
plt.plot(x, lensa.values.T, color='r', alpha = 0.15)
plt.plot(x, lensa.mean('run'), color='k', ls='--', zorder=4)

plt.plot(x, lens.values.T, color='b', alpha = 0.15)
plt.plot(x, lens_ensmean, color='k', zorder=4)

plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([1920, 2100])
plt.legend
plt.savefig('../../../megafires_data/png/LENS_ensemble_anomaly.png', dpi=300)
plt.show()