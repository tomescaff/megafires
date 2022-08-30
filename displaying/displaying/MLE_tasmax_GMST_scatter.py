import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.math as pmath


smf = gmst.get_gmst_annual_5year_smooth()
qnf = stns.get_QN_tmax_jan()

sm = smf.sel(time=slice('1928','2021'))
qn = qnf.sel(time=slice('1928','2021'))

sm = sm.where(sm.time.dt.year != 2017, drop=True)
qn = qn.where(qn.time.dt.year != 2017, drop=True)

xopt = pmath.mle_norm_2d(qn.values, sm.values, [15, 2, 0.5])

mu0, sigma0, alpha = xopt

mu = mu0 + alpha*sm
mu_plus_1sigma = mu + sigma0
mu_plus_2sigma = mu + 2*sigma0

nboot = 10000
filepath = '../../../megafires_data/output/MLE_tasmax_jan_QN_GMST_'+str(nboot)+'.nc'
bspreds = xr.open_dataset(filepath)
bspreds_mu0 = bspreds.mu0.values
bspreds_sigma0 = bspreds.sigma0.values
bspreds_alpha = bspreds.alpha.values

T2017 = smf.sel(time='2017').values
T1880 = smf.sel(time='1880').values

mu_1880 = mu0 + alpha*T1880
mu_2017 = mu0 + alpha*T2017

mu_2017_dist = bspreds_mu0 + bspreds_alpha*T2017
mu_1880_dist = bspreds_mu0 + bspreds_alpha*T1880

mu_2017_inf, mu_2017_sup = np.quantile(mu_2017_dist, [0.025, 0.975], axis = 0)
mu_1880_inf, mu_1880_sup = np.quantile(mu_1880_dist, [0.025, 0.975], axis = 0)

err_mu_2017_inf = mu_2017 - mu_2017_inf
err_mu_2017_sup = mu_2017_sup - mu_2017

err_mu_1880_inf = mu_1880 - mu_1880_inf
err_mu_1880_sup = mu_1880_sup - mu_1880

fig = plt.figure(figsize=(8,5))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10
plt.scatter(sm, qn, s=12, marker='o', edgecolor='blue', facecolor='none')
plt.scatter(smf.sel(time='2017'), qnf.sel(time='2017'), s=12, marker='s', edgecolor='fuchsia', facecolor='none')
plt.plot(sm, mu, color='red', linewidth = 2)
plt.plot(sm, mu_plus_1sigma, color='red', linewidth = 0.5)
plt.plot(sm, mu_plus_2sigma, color='red', linewidth = 0.5)

plt.errorbar(x=T1880, y=mu_1880, yerr=[err_mu_1880_inf, err_mu_1880_sup], lw=1.2, color='r', capsize=3, fmt = '.', capthick=1.5)
plt.errorbar(x=T2017, y=mu_2017, yerr=[err_mu_2017_inf, err_mu_2017_sup], lw=1.2, color='r', capsize=3, fmt = '.', capthick=1.5)

plt.ylim([27, 33.75])
plt.grid(color='grey', lw=0.4, ls='--')
ax = plt.gca()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(direction="in")
plt.xlabel('Global mean surface temperature anomaly (smoothed) [ºC]')
plt.ylabel('January Tmax [ºC]')
plt.savefig('../../../megafires_data/png/MLE_tasmax_GMST_scatter.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()