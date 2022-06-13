import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
import processing.stations as stns

# get time series
qn = ut.get_QN_series().sel(time=slice('1966', '2021'))
cu = ut.get_CU_series().sel(time=slice('1966', '2021'))
ch = ut.get_CH_series().sel(time=slice('1966', '2021'))
regindex = stns.get_regindex_tmax_jan(anom=False).sel(time=slice('1960', '2021'))

# get time series anomalies

qn_anom = qn - qn.sel(time=slice('1981', '2010')).mean('time')
cu_anom = cu - cu.sel(time=slice('1981', '2010')).mean('time')
ch_anom = ch - ch.sel(time=slice('1981', '2010')).mean('time')

# create figure
fig = plt.figure(figsize=(11,6))

# plot the data
plt.plot(cu.time.values, cu.values, lw = 1, marker='o', markersize=3, markeredgecolor='blue', markerfacecolor='white', color='blue', label='Curico')
plt.plot(qn.time.values, qn.values, lw = 1, marker='o', markersize=3, markeredgecolor='r', markerfacecolor='white', color='r', label='Quinta Normal')
plt.plot(ch.time.values, ch.values, lw = 1, marker='o', markersize=3, markeredgecolor='green', markerfacecolor='white', color='k', label='Chillan')
plt.plot(regindex.time.values, regindex.values, lw = 1.5, marker='o', markersize=3, markeredgecolor='k', markerfacecolor='white', color='k', label='RegIndex')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set title and labels
plt.xlabel('Time (years)')
plt.ylabel('January Tmax (ºC)')
plt.savefig('../../../megafires_data/png/stations_time_series_regindex.png', dpi=300)
plt.show()
