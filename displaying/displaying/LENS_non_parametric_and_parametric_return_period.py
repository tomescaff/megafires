import matplotlib.pyplot as plt
# plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
import numpy as np
import xarray as xr
from scipy.stats import norm
from sklearn.utils import resample as bootstrap

# get LENS time series
np_jan_all_runs_1950_2021 = ut.get_LENS_jan_1950_2021_ravel()

# get control run time series
np_jan_control_run = ut.get_control_run_jan()

# fit normal dist
normfit_ar_1950_2021 = norm.fit(np_jan_all_runs_1950_2021)
normfit_cr = norm.fit(np_jan_control_run)

# computing best x ticks
y = np.linspace(20, 34, 1000)
x = 1/norm.sf(y, *normfit_ar_1950_2021)

# computing y values
y_norm_ar_1950_2021 = norm.isf(1/x, *normfit_ar_1950_2021)
y_norm_cr = norm.isf(1/x, *normfit_cr)

#confidence intervals
nboot = 1000

# bootstraping LENS all runs
bspreds = np.zeros((nboot, x.size))
for i in range(nboot):
    z = bootstrap(np_jan_all_runs_1950_2021)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/x, *normfit_)
yinf_ar, ysup_ar = np.quantile(bspreds, [0.025, 0.975], axis = 0)

# bootstraping control run
bspreds = np.zeros((nboot, x.size))
for i in range(nboot):
    z = bootstrap(np_jan_control_run)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/x, *normfit_)
yinf_cr, ysup_cr = np.quantile(bspreds, [0.025, 0.975], axis = 0)


# get non parametric return period
u_LENS, tau_LENS = ut.get_LENS_jan_1950_2021_tau(np_jan_all_runs_1950_2021)
u_CR, tau_CR = ut.get_LENS_jan_1950_2021_tau(np_jan_control_run)

# get anom line
qn = ut.get_QN_series()
qn_2017 = qn.sel(time='2017')
qn_1950_1970_mean = qn.sel(time=slice('1950','1970')).mean('time')

lens = ut.get_LENS_jan()
lens_1950_1970_mean = lens.sel(time=slice('1950','1970')).mean('time')
lens_ev = qn_2017.values - qn_1950_1970_mean.values + lens_1950_1970_mean.values 

# create figure
fig = plt.figure(figsize=(13,7))

# plot the confidence intervals
plt.gca().fill_between(x, yinf_ar, ysup_ar, alpha=.25, label='95% confidence interval ACTUAL')
plt.gca().fill_between(x, yinf_cr, ysup_cr, alpha=.25, label='95% confidence interval COUNTERFACTUAL')

# plot the paramtric curves
plt.plot(x, y_norm_ar_1950_2021, color='b', lw=1.5, alpha = 1, label = 'Parametric return period ACTUAL ', zorder=4)
plt.plot(x, y_norm_cr, color='r', lw=1.5, alpha = 1, ls='--', label = 'Parametric return period CONTRAFACTUAL', zorder=4)

# plot the scatter
plt.scatter(tau_LENS, u_LENS, marker='o', facecolor='lightskyblue', edgecolor='none', color='blue', alpha = 0.5, label = 'Non parametric return period ACTUAL')
plt.scatter(tau_CR, u_CR, marker='o', facecolor='lightcoral', edgecolor='none', color='red', alpha = 0.5, label = 'Non parametric return period CONTRAFACTUAL')

# plot the boxplot with evs
plt.boxplot(lens_ev, notch=False, meanline=True,showmeans=True, positions=[1], patch_artist=True, labels=['1950-1970'], showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))

# plot the ev median value
plt.axhline(np.median(lens_ev), lw=1, color='grey', ls='--', label='EV median value')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=3)

# set x log scale
plt.gca().set_xscale('log')

# set ticks and lims
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
plt.xlim([0.9,11000])
plt.ylim([20, 31])

# set title and labels
plt.xlabel('Return period (years)')
plt.ylabel('January Tmax (ºC)')
plt.title('LENS all runs 1950-2021 (ACTUAL) vs. control run (CONTRAFACTUAL)')
plt.savefig('../../../megafires_data/png/LENS_non_parametric_and_parametric_return_period.png', dpi=300)
plt.show()