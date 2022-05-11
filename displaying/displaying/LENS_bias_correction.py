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

# get LENS jan tmax
lens = da[:, ::12]


# get Quinta Normal time series
qn = ut.get_QN_series()

lens_mean_1950_1970 = lens.sel(time=slice('1950', '1970')).mean('time')
qn_mean_1950_1970 = qn.sel(time=slice('1950', '1970')).mean('time')

lens_anom = lens - lens_mean_1950_1970
lens_corrected = lens_anom + qn_mean_1950_1970

lens_corrected_ensmean = lens_corrected.mean('run')
qmin, qmax = np.quantile(lens_corrected.values, [0.0, 1.0], axis = 0)

x = lens_corrected_ensmean.time.dt.year.values
y = lens_corrected_ensmean.values

fig = plt.figure(figsize=(12,8))
plt.fill_between(x, qmax, qmin, color='grey', alpha = 0.7, label='LENS corrected all runs')
plt.plot(x, y, color='b', label='LENS corrected ensemble mean')
plt.plot(qn.time.dt.year.values, qn.values, color='r', label='Quinta Normal')
plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([1920, 2100])
plt.legend()
plt.savefig('../../../megafires_data/png/LENS_bias_correction_maxmin.png', dpi=300)
plt.show()

fig = plt.figure(figsize=(12,8))
plt.plot(x, lens_corrected.values.T, color='grey', alpha = 0.5)
plt.plot(x, y, color='b', label='LENS corrected ensemble mean')
plt.plot(qn.time.dt.year.values, qn.values, color='r', label='Quinta Normal')
plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([1920, 2100])
plt.legend()
plt.savefig('../../../megafires_data/png/LENS_bias_correction_all_runs.png', dpi=300)
plt.show()

fig = plt.figure(figsize=(12,8))
plt.plot(x, lens_corrected.values.T, color='grey', alpha = 0.5)
plt.plot(x, y, color='b', label='LENS corrected ensemble mean')
plt.plot(qn.time.dt.year.values, qn.values, color='r', label='Quinta Normal')
plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([1950, 2021])
plt.legend()
plt.savefig('../../../megafires_data/png/LENS_bias_correction_all_runs_1950_2021.png', dpi=300)
plt.show()