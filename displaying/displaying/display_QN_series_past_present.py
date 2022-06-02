import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut

# get Quinta Normal time series
da_present = ut.get_QN_series()
da_past = ut.get_QN_past_series()

# create figure
fig = plt.figure(figsize=(12,6))

# # plot the mean value
# plt.axhline(da.mean(), lw=1, color='grey', ls='--', label='mean value')

# plot the data
plt.plot(da_present.time.values, da_present.values, lw = 1, marker='o', markersize=7, markeredgecolor='blue', ls='--', markerfacecolor='white', color='blue', label='Tmax CR2 Explorador')
plt.plot(da_past.time.values, da_past.values, lw = 1, marker='o', markersize=7, markeredgecolor='red', ls='--', markerfacecolor='white', color='red', label='Tmax P. Aceituno')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set title and labels
plt.xlabel('Time (years)')
plt.ylabel('January Tmax (ºC)')
plt.title('January Tmax time series at Quinta Normal')
plt.savefig('../../../megafires_data/png/QN_time_series_past_present.png', dpi=300)
plt.show()
