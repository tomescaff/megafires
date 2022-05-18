import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
import numpy as np
from scipy.stats import norm
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

x_vls = np.arange(24, 30.5, 0.5)
far_vals = np.zeros(x_vls.shape)
rr_vals = np.zeros(x_vls.shape)
far_inf_vals = np.zeros(x_vls.shape)
far_sup_vals = np.zeros(x_vls.shape)
rr_inf_vals = np.zeros(x_vls.shape)
rr_sup_vals = np.zeros(x_vls.shape)


for k, threshold in enumerate(x_vls):

    ##########################################
    # get tau 1
    ##########################################

    # fit normal dist
    normfit = norm.fit(np_jan_all_runs_1950_2021)
    tau_1 = 1/norm.sf(threshold, *normfit)

    ##########################################
    # get tau 0
    ##########################################

    # fit normal dist
    normfit = norm.fit(np_jan_control_run)
    tau_0 = 1/norm.sf(threshold, *normfit)

    rr = tau_0/tau_1
    far = (tau_0-tau_1)/tau_0

    far_vals[k] = far
    rr_vals[k] = rr

    ########################
    # Bootstrap using normal
    ########################

    nboot = 1000#100000#

    rr_norm = np.zeros((nboot,))
    far_norm = np.zeros((nboot,))

    for i in range(nboot):
        ac = bootstrap(np_jan_all_runs_1950_2021)
        cf = bootstrap(np_jan_control_run)
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



# create figure
fig, axs = plt.subplots(2, 1, figsize=(8,7), constrained_layout=True)

plt.sca(axs[0])
plt.gca().fill_between(x_vls, far_inf_vals, far_sup_vals, alpha=.25)
plt.plot(x_vls, far_vals, color='b', lw=1.2)
plt.axvline(lens.sel(time=slice('1950','2021')).mean(), lw=1, ls='--', color='k')
plt.axvline(lens_ev, color='k', lw=1)
plt.grid(lw=0.5, ls='--', color='grey', axis='y')
plt.xlim([24, 30])
plt.ylabel('FAR')

plt.sca(axs[1])
plt.gca().fill_between(x_vls, rr_inf_vals, rr_sup_vals, alpha=.25, color='r')
plt.plot(x_vls, rr_vals, color='r', lw=1.2)
plt.axvline(lens.sel(time=slice('1950','2021')).mean(), lw=1, ls='--', color='k')
plt.axvline(lens_ev, color='k', lw=1)
plt.grid(lw=0.5, ls='--', color='grey', axis='y')
plt.xlim([24, 30])
plt.ylabel('RR')
plt.xlabel('Jan Tmax')

plt.savefig('../../../megafires_data/png/LENS_RR_FAR_anomaly.png', dpi=300)
plt.show()