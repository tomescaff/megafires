from pyexpat.model import XML_CQUANT_REP
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
da_jan = da[:, ::12]
da_jan_ensmean = da_jan.mean('run')
qmin, qmax = np.quantile(da_jan.values, [0.0, 1.0], axis = 0)

basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN_control_run.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan_ctr = da[::12]

x = da_jan_ensmean.time.dt.year.values
y = da_jan_ensmean.values

x_ctr = da_jan_ctr.time.dt.year.values
y_ctr = da_jan_ctr.values

# get Quinta Normal time series
qn = ut.get_QN_series()

fig = plt.figure(figsize=(12,8))

plt.fill_between(x, qmax, qmin, color='grey', alpha = 0.7, label='LENS all runs')
plt.plot(x, y, color='b', label='LENS ensemble mean')
plt.plot(x_ctr, y_ctr, lw=1.2, alpha=1, label='Control run')
plt.plot(qn.time.dt.year.values, qn.values, color='r', label='Quinta Normal')
plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([1920, 2100])
plt.legend()
plt.savefig('../../../megafires_data/png/LENS_tmax_control_run.png', dpi=300)
plt.show()

fig = plt.figure(figsize=(12,8))

plt.fill_between(x, qmax, qmin, color='grey', lw=0.8, alpha = 0.7, label='LENS all runs')
plt.plot(x, y, lw=0.8, color='b', label='LENS ensemble mean')
plt.plot(x_ctr, y_ctr, lw=0.8, alpha=0.5, label='Full control run')
plt.plot(qn.time.dt.year.values, qn.values, lw=0.8, color='r', label='Quinta Normal')
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([400, 2200])
plt.xticks(np.arange(400, 2200, 100), rotation = 0)
plt.legend()
plt.savefig('../../../megafires_data/png/LENS_tmax_full_control_run.png', dpi=300)
plt.show()
