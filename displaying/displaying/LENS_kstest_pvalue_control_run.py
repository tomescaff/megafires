import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys
from scipy import stats

sys.path.append('../../processing')

import processing.utils as ut

basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN_control_run.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan = da[::12]
np_jan_control_run = np.ravel(da_jan.values)

x = np.arange(1920,2101,1)
y = np.zeros(x.shape)


for i, year in enumerate(x):

    basedir = '../../../megafires_data/LENS/'
    filename = 'LENS_tmax_mon_QN.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']-273.15
    da_jan = da[:, ::12]
    da_jan = da_jan.sel(time=slice('1920', str(year)))
    np_jan_all_runs = np.ravel(da_jan.values)


    s, p = stats.ks_2samp(np_jan_all_runs,np_jan_control_run)

    y[i] = p


fig = plt.figure(figsize=(13,8))
plt.plot(x,y, marker='o', markerfacecolor='white', markeredgecolor='blue', color='b', linewidth=1.2)
plt.legend()
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlabel('End year')
plt.ylabel('p-value KS test')
plt.xticks(x[::5], rotation=70)
plt.title(f'P-value of KS test as a function of end year (starting 1920))')
plt.savefig('../../../megafires_data/png/LENS_kstest_pvalue_control_run.png', dpi=300)
plt.show()