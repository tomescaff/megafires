import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
import numpy as np

# get Quinta Normal time series
da = ut.get_QN_series()

# get Quinta Normal return period
u, tau = ut.get_QN_tau()
u2, tau2 = ut.get_QN_tau_remove_max()

# create figure
fig = plt.figure(figsize=(8,6))

# plot the scatter
plt.scatter(tau, u, marker='o', facecolor='lightskyblue', edgecolor='blue', color='blue', alpha = 1, label = 'Non parametric return period')

# plot the scatter
plt.scatter(tau2, u2, marker='o', facecolor='lightcoral', edgecolor='red', color='red', alpha = 1, label = 'Non parametric return period (Jan 2017 removed)')

# plot the 1981-2010 clim
plt.axhline(da.sel(time=slice('1981-01-01','2010-12-31')).mean(), lw=1, color='grey', ls='--', label='1981-2010 mean value')

# plot the 1981-2010 clim
plt.axhline(da.max(), lw=1, color='grey', ls='dotted', label='Jan 2017 value')


# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()
plt.gca().set_xscale('log')
# plt.gca().set_yscale('log')

plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
plt.xlim([0.9,101])
# set title and labels
plt.xlabel('Return period (years)')
plt.ylabel('January Tmax (ºC)')
plt.title('January Tmax return period at Quinta Normal')
plt.savefig('../../../megafires_data/png/QN_non_parametric_return_period.png', dpi=300)
plt.show()
    

