import matplotlib.pyplot as plt
# plt.rcParams["font.family"] = 'Arial'
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys
from sklearn.utils import resample as bootstrap
from scipy.stats import chi2
from scipy.stats import t
sys.path.append('../../processing')

import processing.utils as ut
import processing.lens as lens
import processing.stations as stns


# get Quinta Normal time series
qn = stns.get_QN_tmax_jan()
lens2 = lens.get_LENS2_jan_tmax_QNWE()

qn_1928_2021 = qn.sel(time=slice('1928', '2021'))

lens2_1928_2021 = lens2.sel(time=slice('1928', '2021'))
lens2_std_1928_2021 = lens2_1928_2021.std('time')

########################
# QN std estimation 
########################

std_dev = qn_1928_2021.std(ddof=1)
n = qn_1928_2021.size
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


fig = plt.figure(figsize=(3,6))
plt.boxplot(lens2_std_1928_2021, notch=False, meanline=True,showmeans=True, positions=[0.8], patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.errorbar(x=0.2, y=std_dev, yerr=ci, lw=1.2, color='r', capsize=10, fmt = 'o', capthick=1.5)
plt.grid(ls='--', lw=0.4, color='grey', axis='y')
#plt.yticks(np.arange(26, 31+1.0, 1.0))
plt.ylim([0.3,1.6])
plt.xlim([-0.1,1.1])
plt.xticks([0.2,0.8], ["",""], rotation=0)
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
plt.savefig('../../../megafires_data/png/LENS2_QNWE_bias_std_confinterv.png', dpi=300)
plt.show()
