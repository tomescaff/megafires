import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut

# get Quinta Normal time series
da = ut.get_QN_series()

# get linear trend
dtrend = ut.get_linear_trend()
dtrend_2017r = ut.get_linear_trend_2017r()

# create figure
fig = plt.figure(figsize=(8,6))

# plot the mean value
plt.axhline(da.mean(), lw=1, color='grey', ls='--', label='mean value')

# plot the 1981-2010 clim
plt.axhline(da.sel(time=slice('1981-01-01','2010-12-31')).mean(), lw=1, color='k', ls='dotted', label='1981-2010 mean value')

# plot the data
plt.plot(da.time.values, da.values, lw = 1, marker='o', markersize=7, markeredgecolor='blue', markerfacecolor='white', color='blue', label='Tmax')
plt.plot(da.time.values, da.time.dt.year.values*dtrend['b']+ dtrend['a'], lw = 2.5, alpha=0.7, color='red', label='Linear trend ({:.2f} ºC/dec)'.format(dtrend['b']*10))
plt.plot(da.time.values, da.time.dt.year.values*dtrend_2017r['b']+ dtrend_2017r['a'], lw = 2.5, alpha=1.0, color='orange', label='Linear trend 2017 removed ({:.2f} ºC/dec)'.format(dtrend_2017r['b']*10))

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set title and labels
plt.xlabel('Time (years)')
plt.ylabel('January Tmax (ºC)')
plt.title('January Tmax time series at Quinta Normal')
plt.savefig('../../../megafires_data/png/QN_time_series_with_trend_2017r.png', dpi=300)
plt.show()
