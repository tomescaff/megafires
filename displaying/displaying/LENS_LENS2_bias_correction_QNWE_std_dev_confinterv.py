import matplotlib.pyplot as plt
# plt.rcParams["font.family"] = 'Arial'
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys
from sklearn.utils import resample as bootstrap
from scipy.stats import chi2
sys.path.append('../../processing')

import processing.utils as ut
import processing.lens as lens
import processing.stations as stns

lens1 = lens.get_LENS_jan_tmax_QNWE()
lens2 = lens.get_LENS2_jan_tmax_QNWE()

# get Quinta Normal time series
qn = stns.get_QN_tmax_jan()
qn_mean_1950_1970 = qn.sel(time=slice('1950', '1970')).mean('time')

lens1_mean_1950_1970 = lens1.sel(time=slice('1950', '1970')).mean('time')
lens2_mean_1950_1970 = lens2.sel(time=slice('1950', '1970')).mean('time')

lens1_anom = lens1 - lens1_mean_1950_1970
lens2_anom = lens2 - lens2_mean_1950_1970

lens1_fixed_1950_1970 = lens1_anom + qn_mean_1950_1970
lens2_fixed_1950_1970 = lens2_anom + qn_mean_1950_1970

# 1925-1960
qn_mean_1925_1960 = qn.sel(time=slice('1925', '1960')).mean('time')

lens1_mean_1925_1960 = lens1.sel(time=slice('1925', '1960')).mean('time')
lens2_mean_1925_1960 = lens2.sel(time=slice('1925', '1960')).mean('time')

lens1_anom = lens1 - lens1_mean_1925_1960
lens2_anom = lens2 - lens2_mean_1925_1960

lens1_fixed_1925_1960 = lens1_anom + qn_mean_1925_1960
lens2_fixed_1925_1960 = lens2_anom + qn_mean_1925_1960

#####

init = '1990'
end = '2020'
sig_lens1_1950_1970 = lens1_fixed_1950_1970.sel(time=slice(init, end)).std('time')
sig_lens2_1950_1970 = lens2_fixed_1950_1970.sel(time=slice(init, end)).std('time')
sig_lens1_1925_1960 = lens1_fixed_1925_1960.sel(time=slice(init, end)).std('time')
sig_lens2_1925_1960 = lens2_fixed_1925_1960.sel(time=slice(init, end)).std('time')

qn_seltime = qn.sel(time=slice(init, end))

########################
# QN std estimation 
########################

std_dev = qn_seltime.std(ddof=1)
n = qn_seltime.size
alpha = 0.05
q_l = alpha/2
q_r = 1-alpha/2
chi2_l = chi2.ppf(q_l, n-1)
chi2_r = chi2.ppf(q_r, n-1)
ci_l = np.sqrt((n-1)*std_dev**2/chi2_r)
ci_r = np.sqrt((n-1)*std_dev**2/chi2_l)
ci = np.zeros((2,1))
ci[0] = std_dev-ci_l
ci[1] = ci_r-std_dev

fig = plt.figure(figsize=(11,7))
plt.boxplot(sig_lens1_1950_1970, notch=False, meanline=True,showmeans=True, positions=[1], labels=['LENS1 50-70'], patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(sig_lens1_1925_1960, notch=False, meanline=True,showmeans=True, positions=[2], labels=['LENS1 25-60'],patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(sig_lens2_1950_1970, notch=False, meanline=True,showmeans=True, positions=[3], labels=['LENS2 50-70'],patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(sig_lens2_1925_1960, notch=False, meanline=True,showmeans=True, positions=[4], labels=['LENS2 25-60'],patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.errorbar(x=0, y=std_dev, yerr=ci, lw=1.5, color='b', capsize=10, fmt = 'o', capthick=1.5)

plt.grid(ls='--', lw=0.4, color='grey')
#plt.yticks(np.arange(0.5, 1.5+0.2, 0.2))
plt.ylim([0.3,1.6])
plt.title('Standard deviation for period '+init+'-'+end)
plt.ylabel('Jan Tmax (ºC)')
plt.xlim([-1,5])
plt.xticks(rotation=0)
plt.savefig('../../../megafires_data/png/LENS_bias_correction_std_dev_'+init+'_'+end+'.png', dpi=300)
plt.show()
