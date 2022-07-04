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

# anoms
ev = qn_anom + np.mean(np_jan_all_runs_1991_2020)

tau = 1/norm.sf(ev, *normfit_ar_1991_2020)
past_ev = norm.isf(1/tau, *normfit_ar_1851_1880)
future_ev = norm.isf(1/tau, *normfit_ar_2071_2100)

tau_ee_1851_1880 = 1/norm.sf(ev, *normfit_ar_1851_1880)
tau_ee_2071_2100 = 1/norm.sf(ev, *normfit_ar_2071_2100)

# ########################
# # Bootstrap using normal act vs con
# ########################

# nboot = 100000
# rr_norm = np.zeros((nboot,))
# far_norm = np.zeros((nboot,))

# for i in range(nboot):
#     data_1, data_0 = bootstrap(np_jan_all_runs_1991_2020, np_jan_all_runs_1851_1880)
#     normfit_1 = norm.fit(data_1)
#     normfit_0 = norm.fit(data_0)
#     tau_1_i = 1/norm.sf(ev, *normfit_1)
#     tau_0_i = 1/norm.sf(ev, *normfit_0)
#     rr_norm[i] = tau_0_i/tau_1_i
#     far_norm[i] = (tau_0_i-tau_1_i)/tau_0_i


# tau_1 = 1/norm.sf(ev, *normfit_ar_1991_2020)
# tau_0 = 1/norm.sf(ev, *normfit_ar_1851_1880)
# rr = tau_0/tau_1
# far = (tau_0-tau_1)/tau_0

# yinf_rr, ysup_rr = np.quantile(rr_norm, [0.025, 0.975], axis = 0)
# yinf_far, ysup_far = np.quantile(far_norm, [0.025, 0.975], axis = 0)

# yinf_rr_99 = np.quantile(rr_norm, [0.01], axis = 0)
# yinf_far_99 = np.quantile(far_norm, [0.01], axis = 0)

# print('actual vs. contrafactual')
# print('rr, 95min, 95max, 1inf', rr, yinf_rr, ysup_rr, yinf_rr_99)
# print('far, 95min, 95max, 1inf', far, yinf_far, ysup_far, yinf_far_99)

# act_con_rr = rr_norm
# act_con_far = far_norm

########################
# Bootstrap using normal fut vs act
########################

# nboot = 100000
# rr_norm = np.zeros((nboot,))
# far_norm = np.zeros((nboot,))

# for i in range(nboot):
#     data_1, data_0 = bootstrap(np_jan_all_runs_2071_2100, np_jan_all_runs_1991_2020)
#     normfit_1 = norm.fit(data_1)
#     normfit_0 = norm.fit(data_0)
#     tau_1_i = 1/norm.sf(ev, *normfit_1)
#     tau_0_i = 1/norm.sf(ev, *normfit_0)
#     rr_norm[i] = tau_0_i/tau_1_i
#     far_norm[i] = (tau_0_i-tau_1_i)/tau_0_i


# tau_1 = 1/norm.sf(ev, *normfit_ar_2071_2100)
# tau_0 = 1/norm.sf(ev, *normfit_ar_1991_2020)
# rr = tau_0/tau_1
# far = (tau_0-tau_1)/tau_0

# yinf_rr, ysup_rr = np.quantile(rr_norm, [0.025, 0.975], axis = 0)
# yinf_far, ysup_far = np.quantile(far_norm, [0.025, 0.975], axis = 0)

# yinf_rr_99 = np.quantile(rr_norm, [0.01], axis = 0)
# yinf_far_99 = np.quantile(far_norm, [0.01], axis = 0)

# print('future vs. actual')
# print('rr, 95min, 95max, 1inf', rr, yinf_rr, ysup_rr, yinf_rr_99)
# print('far, 95min, 95max, 1inf', far, yinf_far, ysup_far, yinf_far_99)

# fut_act_rr = rr_norm
# fut_act_far = far_norm

# data1 = np_jan_all_runs_1991_2020
# data0 = np_jan_all_runs_1851_1880
data1 = np_jan_all_runs_2071_2100
data0 = np_jan_all_runs_1991_2020
xmin, xmax = 29, 35
tau_min, tau_max = 10, 10**4
ntaus = 12
x_vls = np.arange(xmin, xmax+0.5, 0.5)
far_vals = np.zeros(x_vls.shape)
rr_vals = np.zeros(x_vls.shape)


far_inf_vals = np.zeros(x_vls.shape)
far_sup_vals = np.zeros(x_vls.shape)
rr_inf_vals = np.zeros(x_vls.shape)
rr_sup_vals = np.zeros(x_vls.shape)

logrange = np.linspace(np.log10(tau_min), np.log10(tau_max), ntaus)
t_vls = 10**logrange
dif_vals = np.zeros(t_vls.shape)
dif_inf_vals = np.zeros(t_vls.shape)
dif_sup_vals = np.zeros(t_vls.shape)

nboot = 100000#1000#

########################
# Bootstrap using normal 
########################

rr_norm_ev = np.zeros((nboot,))
far_norm_ev = np.zeros((nboot,))

for i in range(nboot):
    data_1, data_0 = bootstrap(data1, data0)
    normfit_1 = norm.fit(data_1)
    normfit_0 = norm.fit(data_0)
    tau_1_i = 1/norm.sf(ev, *normfit_1)
    tau_0_i = 1/norm.sf(ev, *normfit_0)
    rr_norm_ev[i] = tau_0_i/tau_1_i
    far_norm_ev[i] = (tau_0_i-tau_1_i)/tau_0_i

normfit_data1 = norm.fit(data1)
normfit_data0 = norm.fit(data0)
tau_1 = 1/norm.sf(ev, *normfit_data1)
tau_0 = 1/norm.sf(ev, *normfit_data0)
rr_ev = tau_0/tau_1
far_ev = (tau_0-tau_1)/tau_0

for k, tau_ in enumerate(t_vls):
    
    ##########################################
    # get temp 1
    ##########################################

    normfit = norm.fit(data1)
    temp_1 = norm.isf(1/tau_, *normfit)

    ##########################################
    # get temp 0
    ##########################################

    normfit = norm.fit(data0)
    temp_0 = norm.isf(1/tau_, *normfit)

    dif_vals[k] = temp_1 - temp_0

    ########################
    # Bootstrap using normal
    ########################

    dif_norm = np.zeros((nboot,))

    for i in range(nboot):
        ac = bootstrap(data1)
        cf = bootstrap(data0)
        normfit_ac = norm.fit(ac)
        normfit_cf = norm.fit(cf)
        temp_1_i = norm.isf(1/tau_, *normfit_ac)
        temp_0_i = norm.isf(1/tau_, *normfit_cf)
        dif_norm[i] = temp_1_i - temp_0_i
    
    dif_inf, dif_sup = np.quantile(dif_norm, [0.025, 0.975], axis = 0)
    dif_inf_vals[k] = dif_inf
    dif_sup_vals[k] = dif_sup

for k, threshold in enumerate(x_vls):

    ##########################################
    # get tau 1
    ##########################################

    # fit normal dist
    normfit = norm.fit(data1)
    tau_1 = 1/norm.sf(threshold, *normfit)

    ##########################################
    # get tau 0
    ##########################################

    # fit normal dist
    normfit = norm.fit(data0)
    tau_0 = 1/norm.sf(threshold, *normfit)

    rr = tau_0/tau_1
    far = (tau_0-tau_1)/tau_0

    far_vals[k] = far
    rr_vals[k] = rr

    ########################
    # Bootstrap using normal
    ########################

    rr_norm = np.zeros((nboot,))
    far_norm = np.zeros((nboot,))

    for i in range(nboot):
        ac = bootstrap(data1)
        cf = bootstrap(data0)
        normfit_ac = norm.fit(ac)
        normfit_cf = norm.fit(cf)
        tau_1_i = 1/norm.sf(threshold, *normfit_ac)
        tau_0_i = 1/norm.sf(threshold, *normfit_cf)
        rr_norm[i] = tau_0_i/tau_1_i
        far_norm[i] = (tau_0_i-tau_1_i)/tau_0_i

    rr_inf, rr_sup = np.quantile(rr_norm, [0.025, 0.975], axis = 0)
    far_inf, far_sup = np.quantile(far_norm, [0.025, 0.975], axis = 0)

    far_inf_vals[k] = far_inf
    far_sup_vals[k] = far_sup
    rr_inf_vals[k] = rr_inf
    rr_sup_vals[k] = rr_sup

# create figure and axes
fig, axs = plt.subplots(1, 4, figsize=(15,3.5), constrained_layout=True)

plt.sca(axs[0])
plt.boxplot(rr_norm_ev, notch=False, meanline=True,showmeans=True, patch_artist=True, whis=(5,95), showfliers=False, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'), labels=[''])
plt.ylabel('Risk Ratio')
rax = axs[0].secondary_yaxis('right', functions = (lambda x: 1-1/x, lambda x: 1/(1-x) ))
rax.set_ylabel('FAR')

ax = axs[0]
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
    tick.label.set_weight('light') 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
    tick.label.set_weight('light')

ax = rax
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
ax.tick_params(axis='both', which='major', labelsize=8)

# plot ax 1
plt.sca(axs[1])
plt.gca().fill_between(t_vls, dif_inf_vals, dif_sup_vals, alpha=.25, color='k')
plt.plot(t_vls, dif_vals, color='k', lw=1.5)
plt.axvline(tau, color='k',  lw=0.8, ls='dotted')
plt.xlim([tau_min, tau_max])
plt.ylim([3.7,4.3])
plt.ylabel('Delta Temperature (ºC)')
plt.xlabel('Return period (yr)')
plt.gca().set_xscale('log')

# plot ax 2
plt.sca(axs[2])
plt.gca().fill_between(x_vls, rr_inf_vals, rr_sup_vals, alpha=.25, color='r')
plt.plot(x_vls, rr_vals, color='r', lw=1.5)
plt.axvline(np.mean(data1), lw=0.8, ls='dotted', color='k')
plt.axvline(ev, color='k',  lw=0.8, ls='dotted')
plt.xlim([xmin, xmax])
plt.ylim([0,250*4])
plt.ylabel('Risk Ratio')
plt.xlabel('January Tmax (ºC)')

# plot ax 3
plt.sca(axs[3])
plt.gca().fill_between(x_vls, far_inf_vals, far_sup_vals, alpha=.25)
plt.plot(x_vls, far_vals, color='b', lw=1.5)
plt.axvline(np.mean(data1), lw=0.8, ls='dotted', color='k')
plt.axvline(ev, color='k',  lw=0.8, ls='dotted')
plt.ylim([0.3,1.0])
plt.xlim([xmin, xmax])
plt.ylabel('FAR')
plt.xlabel('January Tmax (ºC)')

# modify all axes
for ax in axs[1:]:
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.tick_params(direction="in")

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8) 
        tick.label.set_weight('light') 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(8) 
        tick.label.set_weight('light')

plt.savefig('../../../megafires_data/png/LENS2_RR_FAR_future_vs_actual.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()
