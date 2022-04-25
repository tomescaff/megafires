import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
from scipy import signal

# get Quinta Normal time series
da = ut.get_QN_series()

# get linear trend
dtrend = ut.get_linear_trend()

# create figure
fig = plt.figure(figsize=(12,8))

# plot the mean value
plt.axhline(da.mean(), lw=1, color='grey', ls='--', label='mean value')

# plot the 1981-2010 clim
plt.axhline(da.sel(time=slice('1981-01-01','2010-12-31')).mean(), lw=1, color='k', ls='dotted', label='1981-2010 mean value')

# plot the data
plt.plot(da.time.values, da.values, lw = 1, marker='o', markersize=0, markeredgecolor='blue', markerfacecolor='white', color='blue', label='Tmax')
plt.plot(da.time.values, (da.values - dtrend['y_pred']) + da.mean('time').values, marker='o', markersize=0, markeredgecolor='red', markerfacecolor='white', color='red', label='Tmax detrended')
# plt.plot(da.time.values, signal.detrend(da.values) + da.mean('time').values, marker='o', markersize=0, markeredgecolor='orange', markerfacecolor='white', color='orange', label='Tmax detrended')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set title and labels
plt.xlabel('Time (years)')
plt.ylabel('January Tmax (ºC)')
plt.title('January Tmax time series at Quinta Normal')
plt.savefig('../../../megafires_data/png/QN_time_series_linear_detrended.png', dpi=300)
plt.show()
