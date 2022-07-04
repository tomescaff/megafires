import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
import processing.lens as lens
import processing.stations as stns

from scipy.stats import norm
from scipy.stats import genextreme
from scipy.stats import gumbel_r, gumbel_l
from scipy import stats
from scipy.stats import kstest
import numpy as np
import xarray as xr
from sklearn.utils import resample as bootstrap

qn = stns.get_QN_tmax_jan()
qn_1928_2021 = qn.sel(time=slice('1928', '2021'))
qn_mean_1928_2021 = qn_1928_2021.mean('time')
qn_anom = qn.sel(time='2017') - qn.sel(time=slice('1991', '2020')).mean('time')

lens2 = lens.get_LENS2_jan_tmax_QNWE()
lens2_mean_1928_2021 = lens2.sel(time=slice('1928', '2021')).mean('time')
lens2_anom = lens2 - lens2_mean_1928_2021
lens2_fixed_1928_2021 = lens2_anom + qn_mean_1928_2021

# get past period
da_jan_1851_1880 = lens2_fixed_1928_2021.sel(time=slice('1851','1880'))
np_jan_all_runs_1851_1880 = np.ravel(da_jan_1851_1880.values)

# get 1991-2020 period
da_jan_1991_2020 = lens2_fixed_1928_2021.sel(time=slice('1991', '2020'))
np_jan_all_runs_1991_2020 = np.ravel(da_jan_1991_2020.values)

# get 2071-2100 period
da_jan_2071_2100 = lens2_fixed_1928_2021.sel(time=slice('2071', '2100'))
np_jan_all_runs_2071_2100 = np.ravel(da_jan_2071_2100.values)

# fit normal
normfit_ar_1851_1880 = norm.fit(np_jan_all_runs_1851_1880)
normfit_ar_1991_2020 = norm.fit(np_jan_all_runs_1991_2020)
normfit_ar_2071_2100 = norm.fit(np_jan_all_runs_2071_2100)

# computing best x ticks
y = np.linspace(24, 39, 1000)
x = 1/norm.sf(y, *normfit_ar_2071_2100)

# computing y values
y_norm_ar_1851_1880 = norm.isf(1/x, *normfit_ar_1851_1880)
y_norm_ar_1991_2020 = norm.isf(1/x, *normfit_ar_1991_2020)
y_norm_ar_2071_2100 = norm.isf(1/x, *normfit_ar_2071_2100)

#confidence intervals
nboot = 1000

# bootstraping LENS 1851_1880
bspreds = np.zeros((nboot, x.size))
for i in range(nboot):
    z = bootstrap(np_jan_all_runs_1851_1880)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/x, *normfit_)
yinf_1851_1880, ysup_1851_1880 = np.quantile(bspreds, [0.025, 0.975], axis = 0)

# bootstraping LENS 1991-2020
bspreds = np.zeros((nboot, x.size))
for i in range(nboot):
    z = bootstrap(np_jan_all_runs_1991_2020)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/x, *normfit_)
yinf_1991_2020, ysup_1991_2020 = np.quantile(bspreds, [0.025, 0.975], axis = 0)

# bootstraping LENS 2071-2100
bspreds = np.zeros((nboot, x.size))
for i in range(nboot):
    z = bootstrap(np_jan_all_runs_2071_2100)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/x, *normfit_)
yinf_2071_2100, ysup_2071_2100 = np.quantile(bspreds, [0.025, 0.975], axis = 0)

# get non parametric return period
u_1851_1880, tau_1851_1880 = ut.get_LENS_jan_1950_2021_tau(np_jan_all_runs_1851_1880)
u_1991_2020, tau_1991_2020 = ut.get_LENS_jan_1950_2021_tau(np_jan_all_runs_1991_2020)
u_2071_2100, tau_2071_2100 = ut.get_LENS_jan_1950_2021_tau(np_jan_all_runs_2071_2100)

# get anom line

ev = qn_anom + np.mean(np_jan_all_runs_1991_2020)

tau = 1/norm.sf(ev, *normfit_ar_1991_2020)
past_ev = norm.isf(1/tau, *normfit_ar_1851_1880)
future_ev = norm.isf(1/tau, *normfit_ar_2071_2100)

tau_ee_1851_1880 = 1/norm.sf(ev, *normfit_ar_1851_1880)
tau_ee_2071_2100 = 1/norm.sf(ev, *normfit_ar_2071_2100)

# create figure
fig = plt.figure(figsize=(12,7.5))

# plot the confidence intervals
plt.gca().fill_between(x, yinf_1851_1880, ysup_1851_1880, color='b', alpha=.25)
plt.gca().fill_between(x, yinf_1991_2020, ysup_1991_2020, color='k', alpha=.25)
plt.gca().fill_between(x, yinf_2071_2100, ysup_2071_2100, color='r', alpha=.25)

# plot the paramtric curves
plt.plot(x, y_norm_ar_1851_1880, color='b', lw=1.5, alpha = 1, label = 'Parametric return period PAST', zorder=4)
plt.plot(x, y_norm_ar_1991_2020, color='k', lw=1.5, alpha = 1, label = 'Parametric return period PRESENT', zorder=4)
plt.plot(x, y_norm_ar_2071_2100, color='r', lw=1.5, alpha = 1, label = 'Parametric return period FUTURE', zorder=4)

# plot the scatter
#plt.scatter(tau_1851_1880, u_1851_1880, marker='o', facecolor='lightskyblue', edgecolor='none', color='blue', alpha = 0.5)
#plt.scatter(tau_1991_2020, u_1991_2020, marker='o', facecolor='grey', edgecolor='none', color='k', alpha = 0.5)
#plt.scatter(tau_2071_2100, u_2071_2100, marker='o', facecolor='coral', edgecolor='none', color='red', alpha = 0.5)

# plot the ev 
plt.axhline(ev, lw=1, color='grey', ls='-')

# plot the mean values
plt.axhline(np.mean(np_jan_all_runs_1851_1880), lw=1, color='b', ls='dotted')
plt.axhline(np.mean(np_jan_all_runs_1991_2020), lw=1, color='k', ls='dotted')
plt.axhline(np.mean(np_jan_all_runs_2071_2100), lw=1, color='r', ls='dotted')

# plot the ev value
plt.axhline(past_ev, lw=1, color='b', ls='--')
plt.axhline(future_ev, lw=1, color='r', ls='--')

# plot ee tau value
plt.axvline(tau, lw=1, color='k', ls='dotted')
plt.axvline(tau_ee_1851_1880, lw=1, color='b', ls='--')
plt.axvline(tau_ee_2071_2100, lw=1, color='r', ls='--')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set x log scale
plt.gca().set_xscale('log')

# set ticks and lims
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000])
plt.xlim([0.9,10000])
plt.ylim([24.5, 39])

# set title and labels
plt.xlabel('Return period (years)')
plt.ylabel('January Tmax (ºC)')
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
    tick.label.set_weight('light') 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
    tick.label.set_weight('light')

plt.tight_layout()
# plt.savefig('../../../megafires_data/png/LENS2_return_period_past_present_future_normal_full.png', dpi=300)
# plt.show()

# bootstraping LENS 1991-2020
# nboot = 100000
# bspreds = np.zeros((nboot))
# for i in range(nboot):
#     z = bootstrap(np_jan_all_runs_1991_2020)
#     normfit_ = norm.fit(z)
#     bspreds[i] = 1/norm.sf(ev, *normfit_)
# zinf_, zsup_ = np.quantile(bspreds, [0.025, 0.975], axis = 0)

# nboot = 100000
# bspreds = np.zeros((nboot))
# for i in range(nboot):
#     z = bootstrap(np_jan_all_runs_1851_1880)
#     normfit_ = norm.fit(z)
#     bspreds[i] = 1/norm.sf(ev, *normfit_)
# zinf_, zsup_ = np.quantile(bspreds, [0.025, 0.975], axis = 0)

# nboot = 100000
# bspreds = np.zeros((nboot))
# for i in range(nboot):
#     z = bootstrap(np_jan_all_runs_2071_2100)
#     normfit_ = norm.fit(z)
#     bspreds[i] = 1/norm.sf(ev, *normfit_)
# zinf_, zsup_ = np.quantile(bspreds, [0.025, 0.975], axis = 0)