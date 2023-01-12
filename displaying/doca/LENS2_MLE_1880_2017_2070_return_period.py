import sys
import numpy as np
import xarray as xr
import pandas as pd
from scipy.stats import norm

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.lens as lens
import processing.math as pmath

index = ['tau cf', 'tau ac', 'rr c-a', 'far c-a', 'delta c-a']
columns = ['raw', '95ci lower', '95ci upper', '1percentile']
df = pd.DataFrame(columns=columns, index=index)

ac_year = '2017'
cf_year = '1880'

# raw values

lens2_gmst_full = gmst.get_gmst_annual_lens2_ensmean()
lens2_tmax_full = lens.get_LENS2_jan_tmax_QNW()

lens2_gmst = lens2_gmst_full.sel(time=slice('1850', '2021'))
lens2_tmax = lens2_tmax_full.sel(time=slice('1850', '2021'))

lens2_gmst_arr = np.tile(lens2_gmst.values, lens2_tmax.shape[0])
lens2_tmax_arr = np.ravel(lens2_tmax.values)

xopt = pmath.mle_norm_2d(lens2_tmax_arr, lens2_gmst_arr, [29.55, 1.03, 1.11])

mu0, sigma0, alpha = xopt
mu = mu0 + alpha*lens2_gmst_full

mu_MLE_ac = mu.sel(time = ac_year)
mu_MLE_cf = mu.sel(time = cf_year)
sigma_MLE_acf = sigma0

# define tau
tau_ac = 1397

# get ev value 
ev = norm.isf(1/tau_ac, mu_MLE_ac, sigma_MLE_acf)

#################

ac_year = '2017'
fu_year = '2070'

# raw values

lens2_gmst_full = gmst.get_gmst_annual_lens2_ensmean()
lens2_tmax_full = lens.get_LENS2_jan_tmax_QNW()

lens2_gmst = lens2_gmst_full.sel(time=slice('2010', '2100'))
lens2_tmax = lens2_tmax_full.sel(time=slice('2010', '2100'))

lens2_gmst_arr = np.tile(lens2_gmst.values, lens2_tmax.shape[0])
lens2_tmax_arr = np.ravel(lens2_tmax.values)

xopt = pmath.mle_norm_2d(lens2_tmax_arr, lens2_gmst_arr, [29.55, 1.03, 1.11])

mu0, sigma0, alpha = xopt
mu = mu0 + alpha*lens2_gmst_full

mu_MLE_fu = mu.sel(time = fu_year)
sigma_MLE_fu = sigma0

# computing best x ticks
y = np.linspace(24, 39, 1000)
x = 1/norm.sf(y, mu_MLE_fu, sigma_MLE_fu)

# computing y values
y_norm_ar_1851_1880 = norm.isf(1/x, mu_MLE_cf, sigma_MLE_acf)
y_norm_ar_1991_2020 = norm.isf(1/x, mu_MLE_ac, sigma_MLE_acf)
y_norm_ar_2071_2100 = norm.isf(1/x, mu_MLE_fu, sigma_MLE_fu)

# get anom line

ev = norm.isf(1/tau_ac, mu_MLE_ac, sigma_MLE_acf)

tau = 1397
past_ev = norm.isf(1/tau, mu_MLE_cf, sigma_MLE_acf)
future_ev = norm.isf(1/tau, mu_MLE_fu, sigma_MLE_fu)

tau_ee_1851_1880 = 1/norm.sf(ev, mu_MLE_cf, sigma_MLE_acf)
tau_ee_2071_2100 = 1/norm.sf(ev, mu_MLE_fu, sigma_MLE_fu)

# create figure
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
fig = plt.figure(figsize=(12,7.5))

# plot the paramtric curves
plt.plot(x, y_norm_ar_1851_1880, color='b', lw=3.5, alpha = 1, label = 'Parametric return period PAST', zorder=4)
plt.plot(x, y_norm_ar_1991_2020, color='k', lw=3.5, alpha = 1, label = 'Parametric return period PRESENT', zorder=4)
plt.plot(x, y_norm_ar_2071_2100, color='r', lw=3.5, alpha = 1, label = 'Parametric return period FUTURE', zorder=4)

# plot the mean values
plt.axhline(mu_MLE_cf, lw=1, color='b', ls='dotted')
plt.axhline(mu_MLE_ac, lw=1, color='k', ls='dotted')
plt.axhline(mu_MLE_fu, lw=1, color='r', ls='dotted')

# plot the ev value
plt.axhline(ev, lw=1.5, color='k', ls='--')
plt.axhline(past_ev, lw=1.5, color='b', ls='--')
plt.axhline(future_ev, lw=1.5, color='r', ls='--')

# plot ee tau value
plt.axvline(tau, lw=2, color='fuchsia', ls='--')
plt.axvline(tau_ee_1851_1880, lw=2, color='b', ls='--')
plt.axvline(tau_ee_2071_2100, lw=2, color='r', ls='--')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set x log scale
plt.gca().set_xscale('log')

# set ticks and lims
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000])
plt.xlim([0.9,100000])
plt.ylim([24.5, 39])

ax = plt.gca()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(direction="in")


plt.tight_layout()
plt.savefig('../../../megafires_data/png/doca_return_period_present_future_normal.png', dpi=300)
plt.show()