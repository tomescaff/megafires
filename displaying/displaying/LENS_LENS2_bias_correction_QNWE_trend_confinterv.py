import matplotlib.pyplot as plt
# plt.rcParams["font.family"] = 'Arial'
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys
from sklearn.utils import resample as bootstrap
from scipy import stats
from scipy.stats import t
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
lens1_fixed_1950_1970 = lens1_fixed_1950_1970.sel(time=slice(init, end))
lens2_fixed_1950_1970 = lens2_fixed_1950_1970.sel(time=slice(init, end))
lens1_fixed_1925_1960 = lens1_fixed_1925_1960.sel(time=slice(init, end))
lens2_fixed_1925_1960 = lens2_fixed_1925_1960.sel(time=slice(init, end))
qn_seltime = qn.sel(time=slice(init, end))

# compute trend LENS

def get_trends(series):
    nruns, ntime = series.shape
    series_trends = np.zeros((nruns,))
    for k in range(nruns):
        y = series[k,:].values
        x = series[k,:].time.dt.year.values
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        series_trends[k] = slope
    return series_trends

lens1_fixed_1950_1970_trends = get_trends(lens1_fixed_1950_1970)
lens2_fixed_1950_1970_trends = get_trends(lens2_fixed_1950_1970)
lens1_fixed_1925_1960_trends = get_trends(lens1_fixed_1925_1960)
lens2_fixed_1925_1960_trends = get_trends(lens2_fixed_1925_1960)

########################
# compute QN CI trend 
########################

slope, intercept, r_value, p_value, std_err = stats.linregress(qn_seltime.time.dt.year.values, qn_seltime.values)

x = qn_seltime.time.dt.year.values
x_bar = np.mean(x)
n = qn_seltime.size
y = qn_seltime.values
y_hat = x*slope + intercept
SE = np.sqrt(np.sum((y-y_hat)**2)/(n-2))/np.sqrt(np.sum((x-x_bar)**2))

alpha = 0.05
dof = n-2
p_star = 1-alpha/2
ME = t.ppf(p_star, dof)*SE






fig = plt.figure(figsize=(11,7))
plt.boxplot(lens1_fixed_1950_1970_trends*100, notch=False, meanline=True,showmeans=True, positions=[1], labels=['LENS1 50-70'], patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(lens1_fixed_1925_1960_trends*100, notch=False, meanline=True,showmeans=True, positions=[2], labels=['LENS1 25-60'],patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(lens2_fixed_1950_1970_trends*100, notch=False, meanline=True,showmeans=True, positions=[3], labels=['LENS2 50-70'],patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.boxplot(lens2_fixed_1925_1960_trends*100, notch=False, meanline=True,showmeans=True, positions=[4], labels=['LENS2 25-60'],patch_artist=True, showfliers=True, boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
plt.errorbar(x=0, y=slope*100, yerr=ME*100, lw=1.5, color='b', capsize=10, fmt = 'o', capthick=1.5)


plt.grid(ls='--', lw=0.4, color='grey')
plt.yticks(np.arange(-6.5, 10+1, 1.0))
plt.ylim([-6.5,10.0])
plt.title('Trend for period '+init+'-'+end)
plt.ylabel('Jan Tmax trend (ºC/100yr)')
plt.xlim([-1,5])
plt.xticks(rotation=0)
plt.savefig('../../../megafires_data/png/LENS_bias_correction_trend_'+init+'_'+end+'.png', dpi=300)
plt.show()
