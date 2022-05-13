import matplotlib.pyplot as plt
# plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
import numpy as np
import xarray as xr
from scipy.stats import norm

# get LENS time series
basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan = da[:, ::12]
# get 1950-2021 period
da_jan_1950_2021 = da_jan.sel(time=slice('1950', '2021'))
np_jan_all_runs_1950_2021 = np.ravel(da_jan_1950_2021.values)

# get control run time series
basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN_control_run.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan = da[::12]
np_jan_control_run = np.ravel(da_jan.values)

# fit normal dist
normfit_ar_1950_2021 = norm.fit(np_jan_all_runs_1950_2021)
normfit_cr = norm.fit(np_jan_control_run)

x = np.linspace(20, 34, 1000)
y_norm_ar_1950_2021 = 1/norm.sf(x, *normfit_ar_1950_2021)
y_norm_cr = 1/norm.sf(x, *normfit_cr)


# get Quinta Normal return period
# u, tau = ut.get_QN_tau()
# u2, tau2 = ut.get_QN_tau_remove_max()

# create figure
fig = plt.figure(figsize=(12,6))

# plot the paramtric curves
plt.plot(y_norm_ar_1950_2021, x, color='b', lw=1.5, alpha = 1, label = 'Return period LENS all runs 1950-2021 (ACTUAL)')
plt.plot(y_norm_cr, x, color='r', lw=1.5, alpha = 1, label = 'Return period LENS control run (CONTRAFACTUAL)')

# plot the scatter
# plt.scatter(tau, u, marker='o', facecolor='lightskyblue', edgecolor='blue', color='blue', alpha = 1, label = 'Non parametric return period')
# plt.scatter(tau2, u2, marker='o', facecolor='lightcoral', edgecolor='red', color='red', alpha = 1, label = 'Non parametric return period (Jan 2017 removed)')

# plot the 1981-2010 clim
# plt.axhline(da.sel(time=slice('1981-01-01','2010-12-31')).mean(), lw=1, color='grey', ls='--', label='1981-2010 mean value')

# plot the max value line
# plt.axhline(da.max(), lw=1, color='grey', ls='dotted', label='Jan 2017 value')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set x log scale
plt.gca().set_xscale('log')

# set ticks and lims
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
plt.xlim([0.9,11000])
plt.ylim([20, 34])

# set title and labels
plt.xlabel('Return period (years)')
plt.ylabel('January Tmax (ºC)')
plt.title('January Tmax return period LENS all runs vs. control run')
# plt.savefig('../../../megafires_data/png/QN_parametric_return_period.png', dpi=300)
plt.show()