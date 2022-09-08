import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

sys.path.append('../../processing')

import processing.lens as lens
import processing.stations as stns

# get qn data
qn = stns.get_QN_tmax_jan()

# get lens1 data
lens1 = lens.get_LENS_jan_tmax_QNWE()
lens1_sup = lens1.mean('run') + lens1.std('run')
lens1_inf = lens1.mean('run') - lens1.std('run')

# get lens2 data
lens2 = lens.get_LENS2_jan_tmax_QNWE()
lens2_sup = lens2.mean('run') + lens2.std('run')
lens2_inf = lens2.mean('run') - lens2.std('run')

# compute anomaly values
qn_anom = qn - qn.sel(time=slice('1980','2010')).mean('time')
lens1_anom = lens1 - lens1.sel(time=slice('1980','2010')).mean('time')
lens2_anom = lens2 - lens2.sel(time=slice('1980','2010')).mean('time')

lens1_anom_sup = lens1_anom.mean('run') + lens1_anom.std('run')
lens1_anom_inf = lens1_anom.mean('run') - lens1_anom.std('run')

lens2_anom_sup = lens2_anom.mean('run') + lens2_anom.std('run')
lens2_anom_inf = lens2_anom.mean('run') - lens2_anom.std('run')

fig, axs = plt.subplots(3,2, figsize=(15,8))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10

plt.sca(axs[0,0])
plt.plot(lens1.time.dt.year, lens1.values.T, 'skyblue', lw=1)
plt.plot(lens1.time.dt.year, lens1.mean('run'), 'b', label='QN-lens1')
plt.plot(qn.time.dt.year, qn.values, 'r', lw=2, label='QN-obs')
plt.xlim([1840, 2100])
plt.ylim([19, 33])
plt.ylabel('Jan tasmax (ºC)')
plt.legend(loc='upper left')
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.gca().tick_params(direction="in")

plt.sca(axs[1,0])
plt.plot(lens2.time.dt.year, lens2.values.T, 'grey', lw=1)
plt.plot(lens2.time.dt.year, lens2.mean('run'), 'k', label='QN-lens2')
plt.plot(qn.time.dt.year, qn.values, 'r', lw=2, label='QN-obs')
plt.xlim([1840, 2100])
plt.ylim([19, 33])
plt.ylabel('Jan tasmax (ºC)')
plt.legend(loc='upper left')
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.gca().tick_params(direction="in")

plt.sca(axs[2,0])
plt.fill_between(lens1.time.dt.year, lens1_inf, lens1_sup, color='b', lw=1, alpha=0.3, label='QN-lens1 mu+-std')
plt.fill_between(lens2.time.dt.year, lens2_inf, lens2_sup, color='k', lw=1, alpha=0.3, label='QN-lens2 mu+-std')
plt.plot(lens1.time.dt.year, lens1.mean('run'), 'b', label='QN-lens1')
plt.plot(lens2.time.dt.year, lens2.mean('run'), 'k', label='QN-lens2')
plt.plot(qn.time.dt.year, qn.values, 'r', lw=2, label='QN-obs')
plt.xlim([1840, 2100])
plt.ylim([19, 33])
plt.legend(loc='upper left')
plt.ylabel('Jan tasmax (ºC)')
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.gca().tick_params(direction="in")

plt.sca(axs[0,1])
plt.plot(lens1_anom.time.dt.year, lens1_anom.values.T, 'skyblue', lw=1)
plt.plot(lens1_anom.time.dt.year, lens1_anom.mean('run'), 'b')
plt.plot(qn_anom.time.dt.year, qn_anom.values, 'r', lw=2)
plt.xlim([1840, 2100])
plt.ylim([-5, 9])
plt.ylabel('Jan tasmax anom. (ºC)')
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.gca().tick_params(direction="in")

plt.sca(axs[1,1])
plt.plot(lens2_anom.time.dt.year, lens2_anom.values.T, 'grey', lw=1)
plt.plot(lens2_anom.time.dt.year, lens2_anom.mean('run'), 'k')
plt.plot(qn_anom.time.dt.year, qn_anom.values, 'r', lw=2)
plt.xlim([1840, 2100])
plt.ylim([-5, 9])
plt.ylabel('Jan tasmax anom. (ºC)')
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.gca().tick_params(direction="in")

plt.sca(axs[2,1])
plt.fill_between(lens1_anom.time.dt.year, lens1_anom_inf, lens1_anom_sup, color='b', lw=1, alpha=0.3)
plt.fill_between(lens2_anom.time.dt.year, lens2_anom_inf, lens2_anom_sup, color='k', lw=1, alpha=0.3)
plt.plot(lens1_anom.time.dt.year, lens1_anom.mean('run'), 'b')
plt.plot(lens2_anom.time.dt.year, lens2_anom.mean('run'), 'k')
plt.plot(qn_anom.time.dt.year, qn_anom.values, 'r', lw=2)
plt.xlim([1840, 2100])
plt.ylim([-5, 9])
plt.ylabel('Jan tasmax anom. (ºC)')
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.gca().tick_params(direction="in")

plt.tight_layout()
plt.savefig('../../../megafires_data/png/QN_LENS1_LENS2_time_series_validation.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()