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
qn_anom_1952 = qn.sel(time='1952') - qn.sel(time=slice('1931', '1960')).mean('time')
qn_anom_2017 = qn.sel(time='2017') - qn.sel(time=slice('1991', '2020')).mean('time')

lens2 = lens.get_LENS2_jan_tmax_QNWE()
lens2_mean_1928_2021 = lens2.sel(time=slice('1928', '2021')).mean('time')
lens2_anom = lens2 - lens2_mean_1928_2021
lens2_fixed_1928_2021 = lens2_anom + qn_mean_1928_2021

# get 1931-1960 period
da_jan_1931_1960 = lens2_fixed_1928_2021.sel(time=slice('1931','1960')) 
np_jan_all_runs_1931_1960 = np.ravel(da_jan_1931_1960.values) - np.mean(np.ravel(da_jan_1931_1960.values))

# get 1991-2020 period
da_jan_1991_2020 = lens2_fixed_1928_2021.sel(time=slice('1991', '2020'))
np_jan_all_runs_1991_2020 = np.ravel(da_jan_1991_2020.values) - np.mean(np.ravel(da_jan_1991_2020.values))

# fit normal
normfit_ar_1931_1960 = norm.fit(np_jan_all_runs_1931_1960)
normfit_ar_1991_2020 = norm.fit(np_jan_all_runs_1991_2020)

# get taus
tau_1952 = float(1/norm.sf(qn_anom_1952.values, *normfit_ar_1931_1960))
tau_2017 = float(1/norm.sf(qn_anom_2017.values, *normfit_ar_1991_2020))

# get anoms
anom_like_1952_in_1991_2020  = norm.isf(1/tau_1952, *normfit_ar_1991_2020)
anom_like_2017_in_1931_1960  = norm.isf(1/tau_2017, *normfit_ar_1931_1960)

#confidence intervals
nboot = 1000

# bootstraping LENS2 1931-1960 1952 anom
bspreds = np.zeros((nboot,))
for i in range(nboot):
    z = bootstrap(np_jan_all_runs_1931_1960)
    normfit_ = norm.fit(z)
    bspreds[i] = 1/norm.sf(qn_anom_1952.values, *normfit_)
tau_1952_inf, tau_1952_sup = np.quantile(bspreds, [0.025, 0.975], axis = 0)

# bootstraping LENS2 1991-2020 2017 anom
bspreds = np.zeros((nboot,))
for i in range(nboot):
    z = bootstrap(np_jan_all_runs_1991_2020)
    normfit_ = norm.fit(z)
    bspreds[i] = 1/norm.sf(qn_anom_2017.values, *normfit_)
tau_2017_inf, tau_2017_sup = np.quantile(bspreds, [0.025, 0.975], axis = 0)


# create figure
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12,7.5), constrained_layout=True)

xmin, xmax = -5, 5
x = np.linspace(xmin, xmax, 1000)

plt.sca(axs[0])
plt.fill_between(x, 0, norm.pdf(x, *normfit_ar_1931_1960), x >= qn_anom_1952.values, color='grey', label=f'tau={tau_1952:.0f} ({tau_1952_inf:.0f}, {tau_1952_sup:.0f}) yr')
plt.plot(x, norm.pdf(x, *normfit_ar_1931_1960), color='b', label=f'1931-1960 anom dist. (sig={normfit_ar_1931_1960[1]:.2f})')
plt.axvline(qn_anom_1952.values, color='r', label='1952 Tmax anomaly')
plt.axvline(anom_like_2017_in_1931_1960, color='green', label='2017 equivalent Tmax anomaly')


plt.sca(axs[1])
plt.fill_between(x, 0, norm.pdf(x, *normfit_ar_1991_2020), x >= qn_anom_2017.values, color='grey', label=f'tau={tau_2017:.0f} ({tau_2017_inf:.0f}, {tau_2017_sup:.0f}) yr')
plt.plot(x, norm.pdf(x, *normfit_ar_1991_2020), color='b', label=f'1991-2020 anom dist. (sig={normfit_ar_1991_2020[1]:.2f})')
plt.axvline(qn_anom_2017.values, color='r', label='2017 Tmax anomaly')
plt.axvline(anom_like_1952_in_1991_2020, color='green', label='1952 equivalent Tmax anomaly')

for ax in axs:

    plt.sca(ax)
    plt.xlim([xmin, xmax])
    plt.ylim([0, 0.4])
    # set grid
    plt.grid(lw=0.6, ls='--', color='grey')
    plt.legend(loc='upper left')

    # set title and labels
    plt.xlabel('January Tmax (ºC)')
    plt.ylabel('PDF')

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.tick_params(direction="in")

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8) 
        tick.label.set_weight('light') 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(8) 
        tick.label.set_weight('light')

# plt.savefig('../../../megafires_data/png/LENS2_comparison_1952_2017.png', dpi=300)
plt.show()
