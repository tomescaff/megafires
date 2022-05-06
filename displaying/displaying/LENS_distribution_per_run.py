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
da_jan = da_jan.sel(time=slice('1950', '2021'))


qn = ut.get_QN_series()
qn = qn - qn.mean('time')
hist_qn, bins = np.histogram(qn.values, bins=np.arange(-5, 5.25, 0.40), density=True)
width = 1 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

fig = plt.figure(figsize=(13,8))

nruns = 40
matrix = np.zeros((nruns, center.size))

for k in range(nruns):

    da_jan_k = da_jan.sel(run=k)
    da_jan_k_anom = da_jan_k - da_jan_k.mean('time')
    
    # compute histogram
    hist, bins = np.histogram(da_jan_k_anom.values, bins=np.arange(-5, 5.25, 0.40), density=True)
    matrix[k,:] = hist

bp1 = plt.boxplot(matrix, notch=False, meanline=True,showmeans=True, patch_artist=True, widths=0.25, whis=(5,95), showfliers=True, positions=center, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))

sc1 = plt.scatter(center, hist_qn, s=70, edgecolor='k', facecolor='goldenrod', linewidths=0.2, marker='*', alpha = 1, label='QN distribution', zorder=4)

plt.xticks(np.arange(-5, 5.25, 0.40), [x.round(2) for x in np.arange(-5, 5.25, 0.40)], rotation=70)
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlabel('Tmax anomaly')
plt.ylabel('Probability density')
plt.legend([bp1["boxes"][0], sc1], ['LENS distribution (whiskers at p5%-p95%)', 'Quinta Normal distribution'], loc='upper right')
plt.savefig('../../../megafires_data/png/LENS_distribution_per_run.png', dpi=300)
plt.show()
