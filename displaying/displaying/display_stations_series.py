import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut

# get time series
qn = ut.get_QN_series().sel(time=slice('1960', '2021'))
cu = ut.get_CU_series().sel(time=slice('1960', '2021'))
ch = ut.get_CH_series().sel(time=slice('1960', '2021'))

# get time series anomalies

qn_anom = qn - qn.sel(time=slice('1981', '2010')).mean('time')
cu_anom = cu - cu.sel(time=slice('1981', '2010')).mean('time')
ch_anom = ch - ch.sel(time=slice('1981', '2010')).mean('time')

# create figure
fig = plt.figure(figsize=(11,6))

# plot the data
plt.plot(cu_anom.time.values, cu_anom.values, lw = 1, marker='o', markersize=7, markeredgecolor='blue', markerfacecolor='white', color='blue', label='Curico')
plt.plot(qn_anom.time.values, qn_anom.values, lw = 1, marker='o', markersize=7, markeredgecolor='r', markerfacecolor='white', color='r', label='Quinta Normal')
plt.plot(ch_anom.time.values, ch_anom.values, lw = 1, marker='o', markersize=7, markeredgecolor='k', markerfacecolor='white', color='k', label='Chillan')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set title and labels
plt.xlabel('Time (years)')
plt.ylabel('January Tmax (ºC)')
plt.savefig('../../../megafires_data/png/stations_time_series.png', dpi=300)
plt.show()
