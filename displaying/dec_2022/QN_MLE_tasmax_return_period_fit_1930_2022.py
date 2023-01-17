import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.utils as ut
import processing.math as pmath


smf = gmst.get_gmst_annual_5year_smooth_2022()
qnf = stns.get_QN_tmax_dec_1911_2022()

sm = smf.sel(time=slice('1930','2022'))
qn = qnf.sel(time=slice('1930','2022'))

sm = sm.where(sm.time.dt.year != 2022, drop=True)
qn = qn.where(qn.time.dt.year != 2022, drop=True)

xopt = pmath.mle_norm_2d(qn.values, sm.values, [15, 2, 0.5])

mu0, sigma0, alpha = xopt

mu = mu0 + alpha*sm
mu_plus_1sigma = mu + sigma0
mu_plus_2sigma = mu + 2*sigma0

nboot = 10
filepath = '../../../megafires_data/output/MLE_tasmax_jan_QN_GMST_'+str(nboot)+'_normal_validation_dec_2022.nc'
bspreds = xr.open_dataset(filepath)
bspreds_mu0 = bspreds.mu0.values
bspreds_sigma0 = bspreds.sigma0.values
bspreds_alpha = bspreds.alpha.values

T2022 = smf.sel(time='2022').values
T1880 = smf.sel(time='1880').values

mu_1880 = mu0 + alpha*T1880
mu_2022 = mu0 + alpha*T2022

mu_2022_dist = bspreds_mu0 + bspreds_alpha*T2022
mu_1880_dist = bspreds_mu0 + bspreds_alpha*T1880
sigma_dist = bspreds_sigma0

y = np.linspace(24, 34, 1000)
x_norm_2022 = 1/norm.sf(y, mu_2022, sigma0)
x_norm_1880 = 1/norm.sf(y, mu_1880, sigma0)

matrix_2022 = np.zeros((nboot, y.size))
matrix_1880 = np.zeros((nboot, y.size))

for i in range(nboot):
    matrix_2022[i,:] = 1/norm.sf(y, mu_2022_dist[i], sigma_dist[i])
    matrix_1880[i,:] = 1/norm.sf(y, mu_1880_dist[i], sigma_dist[i])

xinf_2022, xsup_2022= np.quantile(matrix_2022, [0.025, 0.975], axis = 0)
xinf_1880, xsup_1880= np.quantile(matrix_1880, [0.025, 0.975], axis = 0)

# get Quinta Normal return period
qn_shift_2022 = qn.values - mu + mu_2022
qn_shift_1880 = qn.values - mu + mu_1880

u_shift_2022, tau_shift_2022 = ut.get_tau(qn_shift_2022.values)
u_shift_1880, tau_shift_1880 = ut.get_tau(qn_shift_1880.values)

fig = plt.figure(figsize=(8,5))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10

# plot the paramtric curves
plt.plot(x_norm_2022, y, color='r', lw=1.2, label = 'Norm shift fit 2022')
plt.plot(x_norm_1880, y, color='b', lw=1.2, label = 'Norm shift fit 1880')

plt.plot(xinf_2022[xinf_2022 > 10], y[xinf_2022 > 10], color='r', lw=0.5)
plt.plot(xsup_2022[xsup_2022 > 10], y[xsup_2022 > 10], color='r', lw=0.5)
plt.plot(xinf_1880[xinf_1880 > 10], y[xinf_1880 > 10], color='b', lw=0.5)
plt.plot(xsup_1880[xsup_1880 > 10], y[xsup_1880 > 10], color='b', lw=0.5)

plt.scatter(tau_shift_2022, u_shift_2022, marker='+', color='r', s=40, lw=0.5)
plt.scatter(tau_shift_1880, u_shift_1880, marker='x', color='b', s=40, lw=0.5)

plt.axhline(qnf.sel(time='2022'), color='fuchsia', lw=1.2, label='Observed 2022')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend(loc='lower right')

# set x log scale
plt.gca().set_xscale('log')

# set ticks and lims
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
plt.xlim([0.9,11000])
plt.ylim([25, 33.75])

# set title and labels
ax=plt.gca()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(direction="in")
plt.xlabel('Return period (yr)')
plt.ylabel('December Tmax (ºC)')
plt.savefig('../../../megafires_data/png/QN_MLE_tasmax_return_period_fit_1930_2022_dec_2022.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()