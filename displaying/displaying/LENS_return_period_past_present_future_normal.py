import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
from scipy.stats import norm
from scipy.stats import genextreme
from scipy.stats import gumbel_r, gumbel_l
from scipy import stats
from scipy.stats import kstest
import numpy as np
import xarray as xr
from sklearn.utils import resample as bootstrap

# get LENS time series
basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan = da[:, ::12]

# get 1920-1950 period
da_jan_1920_1950 = da_jan.sel(time=slice('1920', '1950'))
np_jan_all_runs_1920_1950 = np.ravel(da_jan_1920_1950.values)

# get 1990-2020 period
da_jan_1990_2020 = da_jan.sel(time=slice('1990', '2020'))
np_jan_all_runs_1990_2020 = np.ravel(da_jan_1990_2020.values)

# get 2070-2100 period
da_jan_2070_2100 = da_jan.sel(time=slice('2070', '2100'))
np_jan_all_runs_2070_2100 = np.ravel(da_jan_2070_2100.values)

# fit normal
normfit_ar_1920_1950 = norm.fit(np_jan_all_runs_1920_1950)
normfit_ar_1990_2020 = norm.fit(np_jan_all_runs_1990_2020)
normfit_ar_2070_2100 = norm.fit(np_jan_all_runs_2070_2100)


# computing best x ticks
y = np.linspace(20, 34, 1000)
x = 1/norm.sf(y, *normfit_ar_2070_2100)

# computing y values
y_norm_ar_1920_1950 = norm.isf(1/x, *normfit_ar_1920_1950)
y_norm_ar_1990_2020 = norm.isf(1/x, *normfit_ar_1990_2020)
y_norm_ar_2070_2100 = norm.isf(1/x, *normfit_ar_2070_2100)

#confidence intervals
nboot = 1000

# bootstraping LENS 1920-1950
bspreds = np.zeros((nboot, x.size))
for i in range(nboot):
    z = bootstrap(np_jan_all_runs_1920_1950)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/x, *normfit_)
yinf_1920_1950, ysup_1920_1950 = np.quantile(bspreds, [0.025, 0.975], axis = 0)

# bootstraping LENS 1990-2020
bspreds = np.zeros((nboot, x.size))
for i in range(nboot):
    z = bootstrap(np_jan_all_runs_1990_2020)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/x, *normfit_)
yinf_1990_2020, ysup_1990_2020 = np.quantile(bspreds, [0.025, 0.975], axis = 0)

# bootstraping LENS 2070-2100
bspreds = np.zeros((nboot, x.size))
for i in range(nboot):
    z = bootstrap(np_jan_all_runs_2070_2100)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/x, *normfit_)
yinf_2070_2100, ysup_2070_2100 = np.quantile(bspreds, [0.025, 0.975], axis = 0)

# get non parametric return period
u_1920_1950, tau_1920_1950 = ut.get_LENS_jan_1950_2021_tau(np_jan_all_runs_1920_1950)
u_1990_2020, tau_1990_2020 = ut.get_LENS_jan_1950_2021_tau(np_jan_all_runs_1990_2020)
u_2070_2100, tau_2070_2100 = ut.get_LENS_jan_1950_2021_tau(np_jan_all_runs_2070_2100)

# get anom line
qn = ut.get_QN_series()
qn_2017 = qn.sel(time='2017')
qn_1950_1970_mean = qn.sel(time=slice('1950','1970')).mean('time')

lens = ut.get_LENS_jan()
lens_1950_1970_mean = lens.sel(time=slice('1950','1970')).mean('time')
lens_ev = qn_2017.values - qn_1950_1970_mean.values + lens_1950_1970_mean.values 

tau = 1/norm.sf(np.median(lens_ev), *normfit_ar_1990_2020)
past_ev = norm.isf(1/tau, *normfit_ar_1920_1950)
future_ev = norm.isf(1/tau, *normfit_ar_2070_2100)

tau_ee_1920_1950 = 1/norm.sf(np.median(lens_ev), *normfit_ar_1920_1950)
tau_ee_2070_2100 = 1/norm.sf(np.median(lens_ev), *normfit_ar_2070_2100)

# create figure
fig = plt.figure(figsize=(14,7.5))

# plot the confidence intervals
plt.gca().fill_between(x, yinf_1920_1950, ysup_1920_1950, color='b', alpha=.25)
plt.gca().fill_between(x, yinf_1990_2020, ysup_1990_2020, color='k', alpha=.25)
plt.gca().fill_between(x, yinf_2070_2100, ysup_2070_2100, color='r', alpha=.25)

# plot the paramtric curves
plt.plot(x, y_norm_ar_1920_1950, color='b', lw=1.5, alpha = 1, label = 'Parametric return period PAST', zorder=4)
plt.plot(x, y_norm_ar_1990_2020, color='k', lw=1.5, alpha = 1, label = 'Parametric return period PRESENT', zorder=4)
plt.plot(x, y_norm_ar_2070_2100, color='r', lw=1.5, alpha = 1, label = 'Parametric return period FUTURE', zorder=4)

# plot the scatter
plt.scatter(tau_1920_1950, u_1920_1950, marker='o', facecolor='lightskyblue', edgecolor='none', color='blue', alpha = 0.5)
plt.scatter(tau_1990_2020, u_1990_2020, marker='o', facecolor='grey', edgecolor='none', color='k', alpha = 0.5)
plt.scatter(tau_2070_2100, u_2070_2100, marker='o', facecolor='lightcoral', edgecolor='none', color='red', alpha = 0.5)

# plot the boxplot with evs
plt.boxplot(lens_ev, notch=False, meanline=True,showmeans=True, positions=[1], patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))

# plot the ev median value
plt.axhline(np.median(lens_ev), lw=1, color='grey', ls='--')

# plot the mean values
plt.axhline(np.mean(np_jan_all_runs_1920_1950), lw=1, color='b', ls='-')
plt.axhline(np.mean(np_jan_all_runs_1990_2020), lw=1, color='k', ls='-')
plt.axhline(np.mean(np_jan_all_runs_2070_2100), lw=1, color='r', ls='-')

# plot the ev value
plt.axhline(past_ev, lw=1, color='b', ls='dotted')
plt.axhline(future_ev, lw=1, color='r', ls='dotted')

# plot ee tau value
plt.axvline(tau, lw=1, color='k', ls='dotted')
plt.axvline(tau_ee_1920_1950, lw=1, color='b', ls='dotted')
plt.axvline(tau_ee_2070_2100, lw=1, color='r', ls='dotted')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set x log scale
plt.gca().set_xscale('log')

# set ticks and lims
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000])
plt.xlim([0.9,60000])
plt.ylim([20, 34])

# set title and labels
plt.xlabel('Return period (years)')
plt.ylabel('January Tmax (ºC)')
plt.tight_layout()
plt.savefig('../../../megafires_data/png/LENS_return_period_past_present_future_normal.png', dpi=300)
plt.show()

