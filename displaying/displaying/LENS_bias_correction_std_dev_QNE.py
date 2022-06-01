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

qn_1950_1970 = qn.sel(time=slice('1950', '1970'))
qn_1950_2021 = qn.sel(time=slice('1950', '2021'))

sig_lens_1950_1970 = lens_corrected_1950_1970.std('time')
sig_lens_1950_2021 = lens_corrected_1950_2021.std('time')

sig_qn_1950_1970 = qn_1950_1970.std('time')
sig_qn_1950_2021 = qn_1950_2021.std('time')

fig = plt.figure(figsize=(6,6))
plt.boxplot(sig_lens_1950_1970, notch=False, meanline=True,showmeans=True, positions=[0], patch_artist=True, labels=['1950-1970'], showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(sig_lens_1950_2021, notch=False, meanline=True,showmeans=True, positions=[1], patch_artist=True, labels=['1950-2021'], showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))

plt.scatter([0], sig_qn_1950_1970, marker='*', s=70, edgecolor='k', facecolor='goldenrod', label='Quinta Normal', zorder=4)
plt.scatter([1], sig_qn_1950_2021, marker='*', s=70, edgecolor='k', facecolor='goldenrod', zorder=4)

plt.grid(ls='--', lw=0.4, color='grey')
plt.legend()
#plt.yticks(np.arange(0.5, 1.5+0.2, 0.2))
#plt.ylim([0.5,1.5])
plt.title('LENS standard deviation')
plt.ylabel('Jan Tmax (ºC)')
plt.xlabel('Period')
plt.xlim([-1,2])
plt.savefig('../../../megafires_data/png/LENS_bias_correction_std_dev_QNE.png', dpi=300)
plt.show()
