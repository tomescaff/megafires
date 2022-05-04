import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
import numpy as np
from scipy.stats import norm
from scipy.stats import genextreme
from scipy.stats import gumbel_r
from scipy.optimize import curve_fit
from sklearn.utils import resample as bootstrap

# get Quinta Normal return period
u, tau = ut.get_QN_tau()

# get Quinta Normal time series
da = ut.get_QN_series()

x = np.arange(1.01, 11000, 0.01)

# fit normal dist
normfit = norm.fit(da.values)
y = norm.isf(1/x, *normfit)

#confidence intervals
nboot = 100
bspreds = np.zeros((nboot, x.size))

for i in range(nboot):
    z = bootstrap(da.values)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/x, *normfit_)

yinf, ysup= np.quantile(bspreds, [0.025, 0.975], axis = 0)

# create figure
fig = plt.figure(figsize=(12,6))

# plot the paramtric curves
plt.gca().fill_between(x, ysup, yinf, alpha=.25, label='95% confidence interval')

# plot the scatter
plt.scatter(tau, u, marker='o', facecolor='lightskyblue', edgecolor='blue', color='blue', alpha = 1, label = 'Non parametric return period')

plt.plot(x,y, color='r', lw=0.8, alpha = 1, label = 'Parametric return period (norm)')

# plot the 1981-2010 clim
plt.axhline(da.sel(time=slice('1981-01-01','2010-12-31')).mean(), lw=1, color='grey', ls='--', label='1981-2010 mean value')

# plot the max value line
plt.axhline(da.max(), lw=1, color='grey', ls='dotted', label='Jan 2017 value')

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
plt.title('January Tmax return period at Quinta Normal')
plt.savefig('../../../megafires_data/png/QN_parametric_return_period_normal_bootstrap.png', dpi=300)
plt.show()
    

