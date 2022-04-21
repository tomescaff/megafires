import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut

# get Quinta Normal time series
da = ut.get_QN_series()

# get CR2MET
cr2_near = ut.get_CR2MET_jan().sel(lat=-33.4450, lon=-70.6828, method='nearest')
cr2_mean = ut.get_CR2MET_jan().sel(lat=slice(-37.0, -33.0)).mean(['lat', 'lon'])

# compute anoms
da = da - da.sel(time=slice('1981-01-01','2010-12-31')).mean('time')
cr2_near = cr2_near - cr2_near.sel(time=slice('1981-01-01','2010-12-31')).mean('time') -6
cr2_mean = cr2_mean - cr2_mean.sel(time=slice('1981-01-01','2010-12-31')).mean('time') -12

# create figure
fig = plt.figure(figsize=(8,6))

# plot the mean value
# plt.axhline(da.mean(), lw=1, color='grey', ls='--', label='mean value')

# plot the 1981-2010 clim
# plt.axhline(da.sel(time=slice('1981-01-01','2010-12-31')).mean(), lw=1, color='red', ls='--', label='1981-2010 mean value')

# plot the data
plt.plot(da.time.values, da.values, lw = 1, markersize=7, markeredgecolor='blue', markerfacecolor='white', color='blue', label=' QN observed values\n (-33.45ºS, -70.68ºW)')
plt.plot(cr2_near.time.values, cr2_near.values, lw = 1, markersize=7, markeredgecolor='k', markerfacecolor='white', color='k', label=' CR2MET nearest point\n (-33.42ºS, -70.68ºW)')
plt.plot(cr2_mean.time.values, cr2_mean.values, lw = 1, markersize=7, markeredgecolor='red', markerfacecolor='white', color='red', label=' CR2MET regional average\n (37-33ºS)')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')


# Shrink current axis's height by 10% on the bottom
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# set legend
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3)

# set title and labels
plt.yticks(list(range(-14,4,2)), ['-2','0', '2', '-2', '0', '2', '-2', '0', '2'])
plt.xlabel('Time (years)')
plt.ylabel('January Tmax Anomaly (ºC)')
plt.title('January Tmax Anomaly comparison between Quinta Normal and CR2MET')
plt.savefig('../../../megafires_data/png/QN_CR2MET_time_series_comparison.png', dpi=300)
plt.show()
