import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys

sys.path.append('../../processing')

import processing.utils as ut

basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN.nc'

filepath = basedir + filename

ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15

da_ensmean = da.mean('run')
da_ensmean_annual = da_ensmean.resample(time='1YS').mean('time')
da_annual = da.resample(time='1YS').mean('time')
qmin, qmax = np.quantile(da_annual.values, [0.0, 1.0], axis = 0)

da_jan = da[:, ::12]
da_jan_ensmean = da_jan.mean('run')
qmin, qmax = np.quantile(da_jan.values, [0.0, 1.0], axis = 0)


x = da_jan_ensmean.time.dt.year.values
y = da_jan_ensmean.values

# get Quinta Normal time series
qn = ut.get_QN_series()

fig = plt.figure(figsize=(12,8))

plt.fill_between(x, qmax, qmin, color='grey', alpha = 0.7, label='LENS all runs')
plt.plot(x, y, color='b', label='LENS ensemble mean')
plt.plot(qn.time.dt.year.values, qn.values, color='r', label='Quinta Normal')
plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([1920, 2100])
plt.legend()
plt.savefig('../../../megafires_data/png/LENS_tmax_full_series.png', dpi=300)
plt.show()