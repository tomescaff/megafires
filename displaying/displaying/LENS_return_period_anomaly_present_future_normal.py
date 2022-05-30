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

# get LENS ensemble anomalies
lens_ensmean = da_jan.mean('run')
lens_anom = da_jan - lens_ensmean

# fix LENSA
lens_mean_1920_1940 = da_jan.sel(time=slice('1920', '1940')).mean('time')
lens_anom_mean_1920_1940 = lens_anom.sel(time=slice('1920', '1940')).mean('time')
lensa = lens_anom - lens_anom_mean_1920_1940 + lens_mean_1920_1940

# get 1990-2020 period
da_jan_1990_2020 = da_jan.sel(time=slice('1990', '2020'))
np_jan_all_runs_1990_2020 = np.ravel(da_jan_1990_2020.values)

# get 2070-2100 period
da_jan_2070_2100 = da_jan.sel(time=slice('2070', '2100'))
np_jan_all_runs_2070_2100 = np.ravel(da_jan_2070_2100.values)

# get 1990-2020 period from LENSA
da_lensa_1990_2020 = lensa.sel(time=slice('1990', '2020'))
np_lensa_1990_2020 = np.ravel(da_lensa_1990_2020.values)

# fit normal
normfit_lensa_1990_2020 = norm.fit(np_lensa_1990_2020)
normfit_ar_1990_2020 = norm.fit(np_jan_all_runs_1990_2020)
normfit_ar_2070_2100 = norm.fit(np_jan_all_runs_2070_2100)

# computing best x ticks
y = np.linspace(20, 34, 1000)
x = 1/norm.sf(y, *normfit_ar_2070_2100)

# computing y values
y_norm_lensa_1990_2020 = norm.isf(1/x, *normfit_lensa_1990_2020)
y_norm_ar_1990_2020 = norm.isf(1/x, *normfit_ar_1990_2020)
y_norm_ar_2070_2100 = norm.isf(1/x, *normfit_ar_2070_2100)

#confidence intervals
nboot = 1000

# bootstraping LENSA 1990-2020
bspreds = np.zeros((nboot, x.size))
for i in range(nboot):
    z = bootstrap(np_lensa_1990_2020)
    normfit_ = norm.fit(z)
    bspreds[i] = norm.isf(1/x, *normfit_)
yinf_lensa_1990_2020, ysup_lensa_1990_2020 = np.quantile(bspreds, [0.025, 0.975], axis = 0)

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
u_lensa_1990_2020, tau_lensa_1990_2020 = ut.get_LENS_jan_1950_2021_tau(np_lensa_1990_2020)
u_1990_2020, tau_1990_2020 = ut.get_LENS_jan_1950_2021_tau(np_jan_all_runs_1990_2020)
u_2070_2100, tau_2070_2100 = ut.get_LENS_jan_1950_2021_tau(np_jan_all_runs_2070_2100)

# get anom line
qn = ut.get_QN_series()
anom = qn.sel(time='2017') - qn.sel(time=slice('1990', '2020')).mean('time')
temp_ee = anom.values + np.mean(np_jan_all_runs_1990_2020)

tau = 1/norm.sf(temp_ee, *normfit_ar_1990_2020)
lensa_ev = norm.isf(1/tau, *normfit_lensa_1990_2020)
future_ev = norm.isf(1/tau, *normfit_ar_2070_2100)

tau_ee_lensa_1990_2020 = 1/norm.sf(temp_ee, *normfit_lensa_1990_2020)
tau_ee_2070_2100 = 1/norm.sf(temp_ee, *normfit_ar_2070_2100)

# create figure
fig = plt.figure(figsize=(14,7.5))

# plot the confidence intervals
plt.gca().fill_between(x, yinf_lensa_1990_2020, ysup_lensa_1990_2020, color='b', alpha=.25)
plt.gca().fill_between(x, yinf_1990_2020, ysup_1990_2020, color='k', alpha=.25)
plt.gca().fill_between(x, yinf_2070_2100, ysup_2070_2100, color='r', alpha=.25)

# plot the paramtric curves
plt.plot(x, y_norm_lensa_1990_2020, color='b', lw=1.5, alpha = 1, label = 'Parametric return period PAST', zorder=4)
plt.plot(x, y_norm_ar_1990_2020, color='k', lw=1.5, alpha = 1, label = 'Parametric return period PRESENT', zorder=4)
plt.plot(x, y_norm_ar_2070_2100, color='r', lw=1.5, alpha = 1, label = 'Parametric return period FUTURE', zorder=4)

# plot the scatter
plt.scatter(tau_lensa_1990_2020, u_lensa_1990_2020, marker='o', facecolor='lightskyblue', edgecolor='none', color='blue', alpha = 0.5)
plt.scatter(tau_1990_2020, u_1990_2020, marker='o', facecolor='grey', edgecolor='none', color='k', alpha = 0.5)
plt.scatter(tau_2070_2100, u_2070_2100, marker='o', facecolor='lightcoral', edgecolor='none', color='red', alpha = 0.5)

# plot the ev median value
plt.axhline(temp_ee, lw=1, color='grey', ls='--')

# plot the mean values
plt.axhline(np.mean(np_lensa_1990_2020), lw=1, color='b', ls='-')
plt.axhline(np.mean(np_jan_all_runs_1990_2020), lw=1, color='k', ls='-')
plt.axhline(np.mean(np_jan_all_runs_2070_2100), lw=1, color='r', ls='-')

# plot the ev value
plt.axhline(lensa_ev, lw=1, color='b', ls='dotted')
plt.axhline(future_ev, lw=1, color='r', ls='dotted')

# plot ee tau value
plt.axvline(tau, lw=1, color='k', ls='dotted')
plt.axvline(tau_ee_lensa_1990_2020, lw=1, color='b', ls='dotted')
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
plt.savefig('../../../megafires_data/png/LENS_return_period_anomaly_present_future_normal.png', dpi=300)
plt.show()

