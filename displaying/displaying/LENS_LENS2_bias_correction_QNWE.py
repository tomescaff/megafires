import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys

sys.path.append('../../processing')

import processing.lens as lens
import processing.stations as stns

lens1 = lens.get_LENS_jan_tmax_QNWE()
lens2 = lens.get_LENS2_jan_tmax_QNWE()

# get Quinta Normal time series
qn = stns.get_QN_tmax_jan()
qn_mean_1950_1970 = qn.sel(time=slice('1950', '1970')).mean('time')

lens1_mean_1950_1970 = lens1.sel(time=slice('1950', '1970')).mean('time')
lens2_mean_1950_1970 = lens2.sel(time=slice('1950', '1970')).mean('time')

lens1_anom = lens1 - lens1_mean_1950_1970
lens2_anom = lens2 - lens2_mean_1950_1970

lens1_fixed_1950_1970 = lens1_anom + qn_mean_1950_1970
lens2_fixed_1950_1970 = lens2_anom + qn_mean_1950_1970

lens1_fixed_1950_1970_ensmean = lens1_fixed_1950_1970.mean('run')
lens2_fixed_1950_1970_ensmean = lens2_fixed_1950_1970.mean('run')

# 1925-1960
qn_mean_1925_1960 = qn.sel(time=slice('1925', '1960')).mean('time')

lens1_mean_1925_1960 = lens1.sel(time=slice('1925', '1960')).mean('time')
lens2_mean_1925_1960 = lens2.sel(time=slice('1925', '1960')).mean('time')

lens1_anom = lens1 - lens1_mean_1925_1960
lens2_anom = lens2 - lens2_mean_1925_1960

lens1_fixed_1925_1960 = lens1_anom + qn_mean_1925_1960
lens2_fixed_1925_1960 = lens2_anom + qn_mean_1925_1960

lens1_fixed_1925_1960_ensmean = lens1_fixed_1925_1960.mean('run')
lens2_fixed_1925_1960_ensmean = lens2_fixed_1925_1960.mean('run')



fig, axs = plt.subplots(2, 2, figsize=(14,8))
plt.sca(axs[0,0])
x = lens1_fixed_1950_1970_ensmean.time.dt.year.values
y = lens1_fixed_1950_1970_ensmean.values
plt.plot(x, lens1_fixed_1950_1970.values.T, color='grey', alpha = 0.5)
plt.plot(x, y, color='b', label='LENS1 ensemble mean (fixed 1950-1970')
plt.plot(qn.time.dt.year.values, qn.values, color='r', label='Quinta Normal')
plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([1920, 2100])
plt.xticks(rotation=90)
plt.ylim([25, 38.5])
plt.legend()

plt.sca(axs[0,1])
x = lens2_fixed_1950_1970_ensmean.time.dt.year.values
y = lens2_fixed_1950_1970_ensmean.values
plt.plot(x, lens2_fixed_1950_1970.values.T, color='grey', alpha = 0.5)
plt.plot(x, y, color='b', label='LENS2 ensemble mean (fixed 1950-1970')
plt.plot(qn.time.dt.year.values, qn.values, color='r', label='Quinta Normal')
plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([1920, 2100])
plt.xticks(rotation=90)
plt.ylim([25, 38.5])
plt.legend()

plt.sca(axs[1,0])
x = lens1_fixed_1925_1960_ensmean.time.dt.year.values
y = lens1_fixed_1925_1960_ensmean.values
plt.plot(x, lens1_fixed_1925_1960.values.T, color='grey', alpha = 0.5)
plt.plot(x, y, color='b', label='LENS1 ensemble mean (fixed 1925-1960')
plt.plot(qn.time.dt.year.values, qn.values, color='r', label='Quinta Normal')
plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([1920, 2100])
plt.xticks(rotation=90)
plt.ylim([25, 38.5])
plt.legend()

plt.sca(axs[1,1])
x = lens2_fixed_1925_1960_ensmean.time.dt.year.values
y = lens2_fixed_1925_1960_ensmean.values
plt.plot(x, lens2_fixed_1925_1960.values.T, color='grey', alpha = 0.5)
plt.plot(x, y, color='b', label='LENS2 ensemble mean (fixed 1925-1960')
plt.plot(qn.time.dt.year.values, qn.values, color='r', label='Quinta Normal')
plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
plt.ylabel('Tmax January (ºC)')
plt.grid(ls='--', lw=0.4, color='grey')
plt.xlim([1920, 2100])
plt.xticks(rotation=90)
plt.legend()
plt.ylim([25, 38.5])
plt.tight_layout()
plt.savefig('../../../megafires_data/png/LENS_LENS2_bias_correction_all_runs_QNWE.png', dpi=300)
plt.show()

# fig = plt.figure(figsize=(12,8))
# plt.plot(x, lens_corrected.values.T, color='grey', alpha = 0.5)
# plt.plot(x, y, color='b', label='LENS corrected ensemble mean')
# plt.plot(qn.time.dt.year.values, qn.values, color='r', label='Quinta Normal')
# plt.xticks(np.arange(1920, 2110, 10), rotation = 0)
# plt.ylabel('Tmax January (ºC)')
# plt.grid(ls='--', lw=0.4, color='grey')
# plt.xlim([1950, 2021])
# plt.legend()
# plt.savefig('../../../megafires_data/png/LENS_bias_correction_all_runs_1950_2021_QNE.png', dpi=300)
# plt.show()