import matplotlib.pyplot as plt
# plt.rcParams["font.family"] = 'Arial'
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys

sys.path.append('../../processing')

import processing.utils as ut

basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QNE.nc'

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

lens_corrected_1950_1970 = lens_corrected.sel(time=slice('1950', '1970'))
lens_corrected_1950_2021 = lens_corrected.sel(time=slice('1950', '2021'))
lens_corrected_1990_2020 = lens_corrected.sel(time=slice('1990', '2020'))

qn_1950_1970 = qn.sel(time=slice('1950', '1970'))
qn_1950_2021 = qn.sel(time=slice('1950', '2021'))
qn_1990_2020 = qn.sel(time=slice('1990', '2020'))

mu_lens_1950_1970 = lens_corrected_1950_1970.mean('time')
mu_lens_1950_2021 = lens_corrected_1950_2021.mean('time')
mu_lens_1990_2020 = lens_corrected_1990_2020.mean('time')

mu_qn_1950_1970 = qn_1950_1970.mean('time')
mu_qn_1950_2021 = qn_1950_2021.mean('time')
mu_qn_1990_2020 = qn_1990_2020.mean('time')

fig = plt.figure(figsize=(6,6))
plt.boxplot(mu_lens_1950_1970, notch=False, meanline=True,showmeans=True, positions=[0], patch_artist=True, labels=['1950-1970'], showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(mu_lens_1950_2021, notch=False, meanline=True,showmeans=True, positions=[1], patch_artist=True, labels=['1950-2021'], showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(mu_lens_1990_2020, notch=False, meanline=True,showmeans=True, positions=[2], patch_artist=True, labels=['1990-2020'], showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))

plt.scatter([0], mu_qn_1950_1970, marker='*', s=70, edgecolor='k', facecolor='goldenrod', label='Quinta Normal', zorder=4)
plt.scatter([1], mu_qn_1950_2021, marker='*', s=70, edgecolor='k', facecolor='goldenrod', zorder=4)
plt.scatter([2], mu_qn_1990_2020, marker='*', s=70, edgecolor='k', facecolor='goldenrod', zorder=4)

plt.grid(ls='--', lw=0.4, color='grey')
plt.legend()
plt.yticks(np.arange(29, 31+0.25, 0.25))
plt.ylim([29,31.00])
plt.title('LENS mean value')
plt.ylabel('Jan Tmax (ºC)')
plt.xlabel('Period')
plt.xlim([-1,3])
plt.savefig('../../../megafires_data/png/LENS_bias_correction_mean_QNE.png', dpi=300)
plt.show()
