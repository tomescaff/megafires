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

##########################################
# get Quinta Normal time series original
##########################################

# get Quinta Normal time series
da_orig = ut.get_QN_series()

tmax = da_orig.max('time').values

# fit normal dist
gumfit = gumbel_r.fit(da_orig.values)
tau_1 = 1/gumbel_r.sf(tmax, *gumfit)

##########################################
# get Quinta Normal time series detrended
##########################################

da_det = ut.get_QN_series_detrended()

# fit normal dist
gumfit = gumbel_r.fit(da_det.values)
tau_0 = 1/gumbel_r.sf(tmax, *gumfit)

########################
# Bootstrap using gumbel
########################

nboot = 1000#100000
rr_gum = np.zeros((nboot,))
far_gum = np.zeros((nboot,))

for i in range(nboot):
    orig, det = bootstrap(da_orig.values, da_det)
    gumfit_orig = gumbel_r.fit(orig)
    gumfit_det = gumbel_r.fit(det)
    tau_1_i = 1/gumbel_r.sf(tmax, *gumfit_orig)
    tau_0_i = 1/gumbel_r.sf(tmax, *gumfit_det)
    rr_gum[i] = tau_0_i/tau_1_i
    far_gum[i] = (tau_0_i-tau_1_i)/tau_0_i


########################
# Bootstrap using normal
########################

rr_norm = np.zeros((nboot,))
far_norm = np.zeros((nboot,))

for i in range(nboot):
    orig, det = bootstrap(da_orig.values, da_det)
    normfit_orig = norm.fit(orig)
    normfit_det = norm.fit(det)
    tau_1_i = 1/norm.sf(tmax, *normfit_orig)
    tau_0_i = 1/norm.sf(tmax, *normfit_det)
    rr_norm[i] = tau_0_i/tau_1_i
    far_norm[i] = (tau_0_i-tau_1_i)/tau_0_i



# create figure
fig, axs = plt.subplots(2, 2, figsize=(8,8), constrained_layout=True)

plt.sca(axs[0,0])
plt.boxplot(rr_gum, notch=False, meanline=True,showmeans=True, patch_artist=True, whis=(5,95), showfliers=False, labels=['Gumbel fit'], boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.axhline(1.0, color='k', lw=1)
plt.grid(lw=0.5, ls='--', color='grey', axis='y')
plt.yticks(np.arange(0,5.5,0.5))
plt.ylabel('Risk Ratio')

plt.sca(axs[0,1])
plt.boxplot(rr_norm, notch=False, meanline=True,showmeans=True, patch_artist=True, whis=(5,95), showfliers=False, labels=['Normal fit'], boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.axhline(1.0, color='k', lw=1)
plt.grid(lw=0.5, ls='--', color='grey', axis='y')
plt.yticks(np.arange(0,70,10))
plt.ylabel('Risk Ratio')

plt.sca(axs[1,0])
plt.boxplot(far_gum, notch=False, meanline=True,showmeans=True, patch_artist=True, whis=(5,95), showfliers=False, labels=['Gumbel fit'], boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.axhline(1.0, color='k', lw=1)
plt.grid(lw=0.5, ls='--', color='grey', axis='y')
# plt.yticks(np.arange(0,5.5,0.5))
plt.ylabel('FAR')

plt.sca(axs[1,1])
plt.boxplot(far_norm, notch=False, meanline=True,showmeans=True, patch_artist=True, whis=(5,95), showfliers=False, labels=['Normal fit'], boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.axhline(1.0, color='k', lw=1)
plt.grid(lw=0.5, ls='--', color='grey', axis='y')
# plt.yticks(np.arange(0,5.5,0.5))
plt.ylabel('FAR')

plt.savefig('../../../megafires_data/png/RR_FAR.png', dpi=300)
plt.show()