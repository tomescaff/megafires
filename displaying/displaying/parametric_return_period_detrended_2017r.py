import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
import numpy as np
from scipy.stats import norm
from scipy.stats import genextreme
from scipy.stats import gumbel_r

# get Quinta Normal time series
da = ut.get_QN_series_detrended_2017r()
x = np.linspace(27, 34, 1000)

# fit normal dist
normfit = norm.fit(da.values)
y_norm = 1/norm.sf(x, *normfit)

# fit gev
gevfit = genextreme.fit(da.values)
y_gev = 1/genextreme.sf(x, *gevfit)

# fit gumbel_r
gumrfit = gumbel_r.fit(da.values)
y_gum = 1/gumbel_r.sf(x, *gumrfit)

# get Quinta Normal return period detrended
u, tau = ut.get_QN_tau_detrended_2017r()
# u2, tau2 = ut.get_QN_tau_remove_max()

# create figure
fig = plt.figure(figsize=(12,6))

# plot the paramtric curves
plt.plot(y_norm, x, color='k', lw=0.8, alpha = 1, label = 'Parametric return period (norm)')
plt.plot(y_gev, x, color='g', lw=0.8, alpha = 1, label = 'parametric return period (GEV)')
plt.plot(y_gum, x, color='r', lw=0.8, alpha = 1, label = 'Parametric return period (Gumbel)')


# plot the scatter
plt.scatter(tau, u, marker='o', facecolor='navajowhite', edgecolor='orange', color='orange', alpha = 1, label = 'Non parametric return period')
# plt.scatter(tau2, u2, marker='o', facecolor='lightcoral', edgecolor='red', color='red', alpha = 1, label = 'Non parametric return period (Jan 2017 removed)')

# plot the 1981-2010 clim
plt.axhline(da.sel(time=slice('1981-01-01','2010-12-31')).mean(), lw=1, color='grey', ls='--', label='1981-2010 mean value')

# plot the max value line
plt.axhline(ut.get_QN_series().max(), lw=1, color='grey', ls='dotted', label='Jan 2017 value')
plt.axhline(u.max(), lw=1, color='grey', ls='-.', label='Max value (detrended, 2017 removed)')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set x log scale
plt.gca().set_xscale('log')

# set ticks and lims
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
plt.xlim([0.9,11000])
plt.ylim([27, 34])

# set title and labels
plt.xlabel('Return period (years)')
plt.ylabel('January Tmax (ºC)')
plt.title('January Tmax return period at Quinta Normal (detrended, 2017 removed)')
plt.savefig('../../../megafires_data/png/QN_parametric_return_period_detrended_2017r.png', dpi=300)
plt.show()
    

