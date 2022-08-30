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

fig = plt.figure(figsize=(12,5))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10
x = np.linspace(25,35,100)
plt.fill_between(x, norm.pdf(x, mu_30yr_2017, sigma_30yr_2017), color='red', lw=0.0, alpha=0.1)
plt.fill_between(x, norm.pdf(x, mu_30yr_1943, sigma_30yr_1943), color='blue', lw=0.0, alpha=0.1)

plt.plot(x, norm.pdf(x, mu_MLE_2017, sigma_MLE), color='r', lw=2)
plt.plot(x, norm.pdf(x, mu_MLE_1943, sigma_MLE), color='b', lw=2)
plt.plot(x, norm.pdf(x, mu_MLE_1880, sigma_MLE), color='k', lw=2)

plt.legend(['MLE 2017', 'MLE 1943', 'MLE 1880', '30y 2017', '30y 1943'])
plt.xlabel('January Tmax (ºC)')
plt.ylabel('PDF')
plt.title('Comparing MLE vs. 30yr methodology')
plt.ylim([0, 0.6])
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
plt.savefig('../../../megafires_data/png/comparison_MLE_30yr_dist.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()