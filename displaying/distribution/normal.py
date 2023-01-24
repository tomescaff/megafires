from scipy.stats import norm
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.math as pmath
from scipy.stats import linregress
from scipy.stats import t
np.random.seed(seed=100)

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'Arial'

smf = gmst.get_gmst_annual_5year_smooth_2022().sel(time=slice('1928', '2022'))
qnf = stns.get_QN_tmax_dec_1911_2022().sel(time=slice('1928', '2022'))

mu0, sigma0, alpha = 25.0, 1.0, 1.5

mu = mu0 + alpha*smf


data = np.array([ norm.rvs(loc=float(mu[i]), scale=sigma0, size=1) for i in range(mu.size)])

qnf[:] = data.squeeze()

sm_no2020 = smf.where(smf.time.dt.year != 2020, drop=True)
qn_no2020 = qnf.where(qnf.time.dt.year != 2020, drop=True)

mu1, sigma1 = norm.fit(qn_no2020.values)

xopt = pmath.mle_norm_2d(qn_no2020.values, sm_no2020.values, [25.0, 1.0, 1.5])
mu02, sigma2, alpha2 = xopt
mu2_ = mu02 + alpha2*smf
mu2 = mu2_.sel(time='2020').values

x = np.linspace(15, 35, 500)
fig, axs = plt.subplots(2,1, figsize=(12,6))

plt.sca(axs[0])
plt.fill_between(x, 0, norm.pdf(x, mu1, sigma1), color='lightgrey')
plt.fill_between(x, 0, norm.pdf(x, mu1, sigma1), where= x >= qnf.sel(time='2020').values, color='red')
plt.plot(x, norm.pdf(x, mu1, sigma1), lw=0.5, color='k')
plt.ylim([0,0.5])
plt.xlim([20, 32])
axs[0].spines.right.set_visible(False)
axs[0].spines.top.set_visible(False)
axs[0].tick_params(direction="in")
plt.axvline(qnf.sel(time='2020').values, ls='--', color='fuchsia', lw=1.0)
plt.axvline(mu1, ls='--', color='k', lw=1.0)

plt.sca(axs[1])
plt.fill_between(x, 0, norm.pdf(x, mu2, sigma2), color='lightgrey')
plt.fill_between(x, 0, norm.pdf(x, mu2, sigma2), where= x >= qnf.sel(time='2020').values, color='red')
plt.plot(x, norm.pdf(x, mu2, sigma2), lw=0.5, color='k')
plt.ylim([0,0.5])
plt.xlim([20, 32])
axs[1].spines.right.set_visible(False)
axs[1].spines.top.set_visible(False)
axs[1].tick_params(direction="in")
plt.axvline(qnf.sel(time='2020').values, ls='--', color='fuchsia', lw=1.0)
plt.axvline(mu2, ls='--', color='k', lw=1.0)

plt.show()
