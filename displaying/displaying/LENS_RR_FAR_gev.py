import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
import numpy as np
from scipy.stats import genextreme
from scipy.optimize import curve_fit
from sklearn.utils import resample as bootstrap

# get LENS time series
np_jan_all_runs_1950_2021 = ut.get_LENS_jan_1950_2021_ravel()

# get control run time series
np_jan_control_run = ut.get_control_run_jan()

# get anom line
qn = ut.get_QN_series()
qn_2017 = qn.sel(time='2017')
qn_1950_1970_mean = qn.sel(time=slice('1950','1970')).mean('time')

lens = ut.get_LENS_jan()
lens_1950_1970_mean = lens.sel(time=slice('1950','1970')).mean('time')
lens_ev = np.median(qn_2017.values - qn_1950_1970_mean.values + lens_1950_1970_mean.values)

##########################################
# get tau 1
##########################################

# fit normal dist
gevfit = genextreme.fit(np_jan_all_runs_1950_2021)
tau_1 = 1/genextreme.sf(lens_ev, *gevfit)

##########################################
# get tau 0
##########################################

# fit normal dist
gevfit = genextreme.fit(np_jan_control_run)
tau_0 = 1/genextreme.sf(lens_ev, *gevfit)

# ########################
# # Bootstrap using gumbel
# ########################

nboot = 100#100000##
# rr_gum = np.zeros((nboot,))
# far_gum = np.zeros((nboot,))

# for i in range(nboot):
#     orig, det = bootstrap(da_orig.values, da_det)
#     gumfit_orig = gumbel_r.fit(orig)
#     gumfit_det = gumbel_r.fit(det)
#     tau_1_i = 1/gumbel_r.sf(tmax, *gumfit_orig)
#     tau_0_i = 1/gumbel_r.sf(tmax, *gumfit_det)
#     rr_gum[i] = tau_0_i/tau_1_i
#     far_gum[i] = (tau_0_i-tau_1_i)/tau_0_i

########################
# Bootstrap using normal
########################

rr_gev = np.zeros((nboot,))
far_gev = np.zeros((nboot,))

for i in range(nboot):
    ac = bootstrap(np_jan_all_runs_1950_2021)
    cf = bootstrap(np_jan_control_run)
    gevfit_ac = genextreme.fit(ac)
    gevfit_cf = genextreme.fit(cf)
    tau_1_i = 1/genextreme.sf(lens_ev, *gevfit_ac)
    tau_0_i = 1/genextreme.sf(lens_ev, *gevfit_cf)
    rr_gev[i] = tau_0_i/tau_1_i
    far_gev[i] = (tau_0_i-tau_1_i)/tau_0_i


# create figure
fig, axs = plt.subplots(2, 1, figsize=(4,7), constrained_layout=True)

plt.sca(axs[0])
plt.boxplot(rr_gev, notch=False, meanline=True,showmeans=True, patch_artist=True, whis=(5,95), showfliers=False, labels=['GEV fit'], boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.axhline(1.0, color='k', lw=1)
plt.grid(lw=0.5, ls='--', color='grey', axis='y')
# plt.yticks(np.arange(0,5.5,0.5))
plt.ylabel('Risk Ratio')

# plt.sca(axs[0,1])
# plt.boxplot(rr_norm, notch=False, meanline=True,showmeans=True, patch_artist=True, whis=(5,95), showfliers=False, labels=['Normal fit'], boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
# plt.axhline(1.0, color='k', lw=1)
# plt.grid(lw=0.5, ls='--', color='grey', axis='y')
# plt.yticks(np.arange(0,70,10))
# plt.ylabel('Risk Ratio')

plt.sca(axs[1])
plt.boxplot(far_gev, notch=False, meanline=True,showmeans=True, patch_artist=True, whis=(5,95), showfliers=False, labels=['GEV fit'], boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.axhline(1.0, color='k', lw=1)
plt.grid(lw=0.5, ls='--', color='grey', axis='y')
# plt.yticks(np.arange(0,5.5,0.5))
plt.ylabel('FAR')

# plt.sca(axs[1,1])
# plt.boxplot(far_norm, notch=False, meanline=True,showmeans=True, patch_artist=True, whis=(5,95), showfliers=False, labels=['Normal fit'], boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
# plt.axhline(1.0, color='k', lw=1)
# plt.grid(lw=0.5, ls='--', color='grey', axis='y')
# # plt.yticks(np.arange(0,5.5,0.5))
# plt.ylabel('FAR')

plt.savefig('../../../megafires_data/png/LENS_RR_FAR_gev.png', dpi=300)
plt.show()