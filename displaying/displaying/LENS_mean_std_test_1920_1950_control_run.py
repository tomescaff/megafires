import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
from scipy.stats import norm
from scipy.stats import genextreme
from scipy.stats import gumbel_r, gumbel_l
from scipy import stats
from scipy.stats import kstest
import numpy as np
import xarray as xr

# get LENS time series
basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan_lens = da[:, ::12]

# get 1920-1950 period
da_jan_1920_1950 = da_jan_lens.sel(time=slice('1920', '1950'))
np_jan_all_runs_1920_1950 = np.ravel(da_jan_1920_1950.values)

# get control run time series
basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN_control_run.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan_cr = da[::12]
np_jan_control_run = np.ravel(da_jan_cr.values)

t_mean, p_mean = stats.ttest_ind(np_jan_all_runs_1920_1950, np_jan_control_run, equal_var=False)
t_var, p_var = stats.levene(np_jan_all_runs_1920_1950, np_jan_control_run)

x = np.arange(1920, 2100, 1)
mat_p_mean = np.zeros(x.shape)
mat_p_var = np.zeros(x.shape)

for i, year in enumerate(x):
    lens_i = np.ravel(da_jan_lens.sel(time=slice('1920', str(year))).values)
    t_mean_i, p_mean_i = stats.ttest_ind(lens_i, np_jan_control_run, equal_var=False)
    t_var_i, p_var_i = stats.levene(lens_i, np_jan_control_run)
    mat_p_mean[i] = p_mean_i
    mat_p_var[i] = p_var_i

fig = plt.figure(figsize=(12,6))
plt.plot(x, mat_p_mean, color='b', lw=1.2, label='p-value mean test')
plt.plot(x, mat_p_var, color='r', lw=1.2, label='p-value var test')
plt.axhline(0.1, color='k', ls='--')
plt.ylim([0,1])
plt.xlim([1920, 2100])
plt.xlabel('End year')
plt.ylabel('p-value')
plt.grid(ls='--', color='gray', lw=0.7)
plt.legend()
plt.savefig('../../../megafires_data/png/LENS_mean_std_test.png', dpi=300)
plt.show()


