import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys
from scipy import stats

sys.path.append('../../processing')

import processing.utils as ut

basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN.nc'

filepath = basedir + filename

ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan = da[:, ::12]
da_jan = da_jan.sel(time=slice('1950', '2021'))
da_jan = da_jan - da_jan.mean('time')
np_jan_all_runs = np.ravel(da_jan.values)
hist_lens, bins = np.histogram(np_jan_all_runs, bins=np.arange(-5, 5.25, 0.4), density=True)

qn = ut.get_QN_series()
qn = qn - qn.mean('time')
hist_qn, bins = np.histogram(qn.values, bins=np.arange(-5, 5.25, 0.4), density=True)

width = 1 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

s, p = stats.ks_2samp(np_jan_all_runs,qn.values)

fig = plt.figure(figsize=(13,8))
plt.bar(center, hist_lens, align='center', width=width, edgecolor='k', facecolor='b', color='blue', alpha = 0.25, label='LENS')
plt.bar(center, hist_qn, align='center', width=width, edgecolor='k', facecolor='r', color='blue', alpha = 0.25, label='QN')
plt.legend()
plt.grid(ls='--', lw=0.4, color='grey')
plt.xticks(np.arange(-5, 5.25, 0.40), [x.round(2) for x in np.arange(-5, 5.25, 0.40)], rotation=70)
plt.xlabel('Tmax anomaly')
plt.ylabel('Probability density')
plt.title(f'LENS (all runs) vs. QN distributions (KS test: p={p:.2f})')
plt.savefig('../../../megafires_data/png/LENS_distribution_all_runs.png', dpi=300)
plt.show()