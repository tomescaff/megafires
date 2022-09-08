import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import t
from scipy.stats import chi2
from scipy.stats import linregress

sys.path.append('../../processing')

import processing.lens as lens
import processing.stations as stns

# get qn data
qn = stns.get_QN_tmax_jan()

# get lens1 data
lens1 = lens.get_LENS_jan_tmax_QNWE()

# get lens2 data
lens2 = lens.get_LENS2_jan_tmax_QNWE()

# common period
qn = qn.sel(time=slice('1930','2021'))
lens1 = lens1.sel(time=slice('1930','2021'))
lens2 = lens2.sel(time=slice('1930','2021'))

# qn mu CI
mean = qn.mean()
std_dev = qn.std(ddof=1)
n = qn.size
dof = n-1
alpha = 0.05
p_star = 1-alpha/2
t_star = t.ppf(p_star, dof)
ME = t_star*std_dev/np.sqrt(n)

# qn std CI
std_dev = qn.std(ddof=1)
n = qn.size
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

# qn trend CI
x = qn.time.dt.year.values
y = qn.values
slope, intercept, r_value, p_value, std_err = linregress(x, y)
y_hat = x*slope + intercept
x_bar = np.mean(x)
n = qn.size
SE = np.sqrt(np.sum((y-y_hat)**2)/(n-2))/np.sqrt(np.sum((x-x_bar)**2))
alpha = 0.05
dof = n-2
p_star = 1-alpha/2
ME_trend = t.ppf(p_star, dof)*SE

def get_trends(series):
    nruns, ntime = series.shape
    series_trends = np.zeros((nruns,))
    for k in range(nruns):
        y = series[k,:].values
        x = series[k,:].time.dt.year.values
        series_trends[k] = linregress(x,y).slope
    return series_trends

lens1_trends = get_trends(lens1)
lens2_trends = get_trends(lens2)

# plot
fig, axs = plt.subplots(1,3, figsize=(15,7))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10

plt.sca(axs[0])
plt.errorbar(x=0.2, y=mean, yerr=ME, lw=1.2, color='r', capsize=10, fmt = 'o', capthick=1.5)
plt.boxplot(lens1.mean('time'), notch=False, meanline=True,showmeans=True, positions=[0.8], patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(lens2.mean('time'), notch=False, meanline=True,showmeans=True, positions=[1.4], patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.ylim([22.5, 30.2])
plt.xlim([-0.1,1.7])
plt.ylabel('January tasmax at QN - 1930-2021 mean value (ºC)')
plt.xticks([0.2,0.8, 1.4], ["DMC (95% CI)","LENS1", "LENS2"], rotation=0)
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.gca().tick_params(direction="in")
plt.title('Mean value')

plt.sca(axs[1])
plt.errorbar(x=0.2, y=std_dev, yerr=ci, lw=1.2, color='r', capsize=10, fmt = 'o', capthick=1.5)
plt.boxplot(lens1.std('time'), notch=False, meanline=True,showmeans=True, positions=[0.8], patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(lens2.std('time'), notch=False, meanline=True,showmeans=True, positions=[1.4], patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.ylim([0.7,1.5])
plt.xlim([-0.1,1.7])
plt.ylabel('January tasmax at QN - 1930-2021 standard deviation (ºC)')
plt.xticks([0.2,0.8, 1.4], ["DMC (95% CI)","LENS1", "LENS2"], rotation=0)
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.gca().tick_params(direction="in")
plt.title('Standard deviation')

plt.sca(axs[2])
plt.errorbar(x=0.2, y=slope*100, yerr=ME_trend*100, lw=1.2, color='r', capsize=10, fmt = 'o', capthick=1.5)
plt.boxplot(lens1_trends*100, notch=False, meanline=True,showmeans=True, positions=[0.8], patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(lens2_trends*100, notch=False, meanline=True,showmeans=True, positions=[1.4], patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.ylim([-1,4])
plt.xlim([-0.1,1.7])
plt.ylabel('January tasmax at QN - 1930-2021 linear trend (ºC/100yr)')
plt.xticks([0.2,0.8, 1.4], ["DMC (95% CI)","LENS1", "LENS2"], rotation=0)
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.gca().tick_params(direction="in")
plt.title('Linear trend')

plt.tight_layout()
plt.savefig('../../../megafires_data/png/QN_LENS1_LENS2_stats_validation.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()