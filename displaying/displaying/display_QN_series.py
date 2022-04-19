import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut

# get Quinta Normal time series
da = ut.get_QN_series()

# create figure
fig = plt.figure(figsize=(8,6))

# plot the mean value
plt.axhline(da.mean(), lw=1, color='grey', ls='--', label='mean value')

# plot the data
plt.plot(da.time.values, da.values, lw = 1, marker='o', markersize=7, markeredgecolor='blue', markerfacecolor='white', color='blue', label='Tmax')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set title and labels
plt.xlabel('Time (years)')
plt.ylabel('January Tmax (ºC)')
plt.title('January Tmax time series at Quinta Normal')
plt.savefig('../../../megafires_data/png/QN_time_series.png', dpi=300)
plt.show()
