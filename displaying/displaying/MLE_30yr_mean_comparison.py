import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.math as pmath

qn = stns.get_QN_tmax_jan().sel(time=slice('1928','2021'))

smf = gmst.get_gmst_annual_5year_smooth()
qnf = stns.get_QN_tmax_jan()

sm = smf.sel(time=slice('1928','2021'))
qn = qnf.sel(time=slice('1928','2021'))

sm_no2017 = sm.where(sm.time.dt.year != 2017, drop=True)
qn_n02017 = qn.where(qn.time.dt.year != 2017, drop=True)

xopt = pmath.mle_norm_2d(qn_n02017.values, sm_no2017.values, [15, 2, 0.5])

mu0, sigma0, alpha = xopt
mu = mu0 + alpha*smf

mu_MLE_2017 = mu.sel(time='2017')
mu_MLE_1943 = mu.sel(time='1943')
mu_MLE_1880 = mu.sel(time='1880')
sigma_MLE = sigma0

qn_30yr_2017 = qnf.sel(time=slice('1992','2021'))
qn_30yr_1943 = qnf.sel(time=slice('1928','1957'))
mu_30yr_2017, sigma_30yr_2017 = norm.fit(qn_30yr_2017)
mu_30yr_1943, sigma_30yr_1943 = norm.fit(qn_30yr_1943)

# bootstrap MLE
nboot = 10000
filepath = '../../../megafires_data/output/MLE_tasmax_jan_QN_GMST_'+str(nboot)+'.nc'
bspreds = xr.open_dataset(filepath)
bspreds_mu0 = bspreds.mu0.values
bspreds_sigma0 = bspreds.sigma0.values
bspreds_alpha = bspreds.alpha.values

T2017 = smf.sel(time='2017').values
T1943 = smf.sel(time='1943').values
T1880 = smf.sel(time='1880').values

mu_2017_dist = bspreds_mu0 + bspreds_alpha*T2017
mu_1943_dist = bspreds_mu0 + bspreds_alpha*T1943
mu_1880_dist = bspreds_mu0 + bspreds_alpha*T1880

mu_mle_2017_inf, mu_mle_2017_sup = np.quantile(mu_2017_dist, [0.025, 0.975], axis = 0)
mu_mle_1943_inf, mu_mle_1943_sup = np.quantile(mu_1943_dist, [0.025, 0.975], axis = 0)
mu_mle_1880_inf, mu_mle_1880_sup = np.quantile(mu_1880_dist, [0.025, 0.975], axis = 0)

# bootstrap 30yr
bspreds_mu_2017 = np.zeros((nboot,))
bspreds_mu_1943 = np.zeros((nboot,))
for i in range(nboot):
    z2017_i = bootstrap(qnf.sel(time=slice('1992','2021')).values)
    z1943_i = bootstrap(qnf.sel(time=slice('1928','1957')).values)
    mu_2017_i, sigma_2017_i = norm.fit(z2017_i) 
    mu_1943_i, sigma_1943_i = norm.fit(z1943_i)
    bspreds_mu_2017[i] = mu_2017_i
    bspreds_mu_1943[i] = mu_1943_i

mu_30y_2017_inf, mu_30y_2017_sup = np.quantile(bspreds_mu_2017, [0.025, 0.975], axis = 0)
mu_30y_1943_inf, mu_30y_1943_sup = np.quantile(bspreds_mu_1943, [0.025, 0.975], axis = 0)

fig = plt.figure(figsize=(12,6))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10
plt.plot(qn.time.dt.year, qn.values, color='grey', lw=0.8, marker='.', ls='--')
plt.plot(mu.time.dt.year, mu.values, color='r', lw=0.8)
qn_roll = qn.rolling(time=30, min_periods=30, center=True).mean('time').dropna('time')
plt.plot(qn_roll.time.dt.year, qn_roll.values, color='b', lw=0.8)
plt.errorbar(x=1943+0.2, y=mu_MLE_1943.values, yerr=[mu_MLE_1943-mu_mle_1943_inf, mu_mle_1943_sup-mu_MLE_1943], lw=1.2, color='r', capsize=6, fmt = '.', capthick=1.5)
plt.errorbar(x=1943-0.2, y=mu_30yr_1943, yerr=[[mu_30yr_1943-mu_30y_1943_inf], [mu_30y_1943_sup-mu_30yr_1943]], lw=1.2, color='b', capsize=6, fmt = '.', capthick=1.5)
plt.errorbar(x=2017+0.2, y=mu_MLE_2017.values, yerr=[mu_MLE_2017-mu_mle_2017_inf, mu_mle_2017_sup-mu_MLE_2017], lw=1.2, color='r', capsize=6, fmt = '.', capthick=1.5)
plt.errorbar(x=2017-0.2, y=mu_30yr_2017, yerr=[[mu_30yr_2017-mu_30y_2017_inf], [mu_30y_2017_sup-mu_30yr_2017]], lw=1.2, color='b', capsize=6, fmt = '.', capthick=1.5)
plt.errorbar(x=1880+0.2, y=mu_MLE_1880.values, yerr=[mu_MLE_1880-mu_mle_1880_inf, mu_mle_1880_sup-mu_MLE_1880], lw=1.2, color='r', capsize=6, fmt = '.', capthick=1.5)
plt.xlabel('Time')
plt.ylabel('January Tmax (ºC)')
plt.legend(['QN Jan Tmax', 'mu=mu(t) -- MLE', 'mu=mu(t) -- 30yr rolling average', 'mu 95%CI -- MLE', 'mu 95%CI -- 30yr'], frameon=False)
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
plt.savefig('../../../megafires_data/png/comparison_MLE_30yr_mean.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()
