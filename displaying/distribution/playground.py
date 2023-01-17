from scipy.stats import norm
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.math as pmath
from scipy.stats import linregress
from scipy.stats import t
np.random.seed(seed=100)
smf = gmst.get_gmst_annual_5year_smooth_2022().sel(time=slice('1928', '2022'))
qnf = stns.get_QN_tmax_dec_1911_2022().sel(time=slice('1928', '2022'))

mu0, sigma0, alpha = 25.0, 1.0, 1.5

mu = mu0 + alpha*smf


data = np.array([ norm.rvs(loc=float(mu[i]), scale=sigma0, size=1) for i in range(mu.size)])

qnf[:] = data.squeeze()
#qnf[-1] = 28.5
sm = smf
qn = qnf

# get linear trend
x = qn.time.dt.year.values
y = qn.values
slope, intercept, rvalue, pvalue, stderr  = linregress(x,y)
trend = x*slope + intercept

# compute confidence intervals
xmean = np.mean(x)
SXX = np.sum((x-xmean)**2)
SSE = np.sum((y - trend)**2)
n = x.size
sig_hat_2_E = SSE/(n-2)
sig_hat_E = np.sqrt(sig_hat_2_E)
rad = np.sqrt(1/n + (x-xmean)**2/SXX) 
p = 0.95
q = (1+p)/2
dof = n-2
tval = t.ppf(q, dof)
delta = tval*sig_hat_E*rad
y_sup = trend+delta
y_inf = trend-delta

# create figure
fig = plt.figure(figsize=(15,6))
gs = fig.add_gridspec(1, 3,  width_ratios=(1, 7, 1),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.2, hspace=0.07)
ax_line = fig.add_subplot(gs[0, 1])
ax_norm = fig.add_subplot(gs[0, 0], sharey=ax_line)
ax_fet = fig.add_subplot(gs[0, 2], sharey=ax_line)

# plot the mean value
plt.sca(ax_line)
plt.plot(x, x*0+qn.mean().values, lw=1.5, color='b', ls='--')
plt.plot(qn.sel(time=slice('1928','1957')).time.dt.year, qn.sel(time=slice('1928','1957')).time.dt.year*0+qn.sel(time=slice('1928','1957')).mean().values, lw=1.5, color='green', ls='--')
plt.plot(qn.sel(time=slice('1992','2021')).time.dt.year, qn.sel(time=slice('1992','2021')).time.dt.year*0+qn.sel(time=slice('1992','2021')).mean().values, lw=1.5, color='brown', ls='--')

print('first 30 years mean', qn.sel(time=slice('1928','1957')).mean())
print('last 30 years mean', qn.sel(time=slice('1992','2021')).mean())

# plot the data
plt.plot(x, y, lw = 1, marker='o', markersize=7, markeredgecolor='blue', markerfacecolor='white', color='blue')
plt.plot(x, trend, lw = 1.5, alpha=0.7, color='r') 
plt.fill_between(x, y_sup, y_inf, facecolor='grey', linewidth=2, alpha=0.2)
plt.plot(qn.sel(time='2020').time.dt.year, qn.sel(time='2020').values, lw=1, marker='o', markersize=7, markeredgecolor='blue', markerfacecolor='blue', color='blue', alpha=0.4)

print('Linear trend ºC/dec)', slope*10)
print('pvalue)', pvalue)

# set grid
plt.grid(lw=0.4, ls='--', color='grey')

# set title and labels
plt.xlabel('Time (year)')

# Hide the right and top spines
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(9) 
    tick.label.set_weight('light') 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(9) 
    tick.label.set_weight('light')

past_climate = qn.sel(time=slice('1928','1957'))
psnt_climate = qn.sel(time=slice('1992','2021'))
psnt_climate_year = psnt_climate.time.dt.year
psnt_climate_no_2017 = psnt_climate.where(psnt_climate_year !=2017, drop=True)
normfit_28_57 = norm.fit(past_climate.values)
# normfit_92_21 = norm.fit(psnt_climate_no_2017.values)
normfit_92_21 = norm.fit(psnt_climate.values)

plt.sca(ax_norm)
ax = plt.gca()
ymin, ymax = ax.get_ylim()
xval = np.linspace(ymin-0.5, ymax, 100)
plt.fill_betweenx(  xval, xval*0, norm.pdf(xval, *normfit_28_57), facecolor='green', linewidth=2, alpha=0.2)
plt.xlim([0, norm.pdf(xval, *normfit_28_57).max()])
plt.fill_betweenx(  xval, xval*0, norm.pdf(xval, *normfit_92_21), facecolor='brown', linewidth=2, alpha=0.2)
plt.xlim([0, norm.pdf(xval, *normfit_92_21).max()])
plt.axhline(qn.sel(time=slice('1928','1957')).mean(), color='green', lw=1.5, ls='--')
plt.axhline(qn.sel(time=slice('1992','2021')).mean(), color='brown', lw=1.5, ls='--')
plt.axhline(qn.sel(time='2020'), color='b', lw=1.5)#, ls='--')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(0) 
    tick.label.set_color('white') 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(9) 
    tick.label.set_weight('light')
plt.ylabel('Tmax December (ºC)')

plt.sca(ax_fet)
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(0) 
    tick.label.set_color('white') 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(9) 
    tick.label.set_weight('light')
plt.bar([0], height=2020*slope + intercept-qn.mean().values, bottom=qn.mean().values, width=0.75, color='r', alpha=0.4)
plt.bar([1.35], height=qn.sel(time='2020').values-qn.mean().values, bottom=qn.mean().values, width=0.75, color='b', alpha=0.4)
plt.axhline(qn.mean().values,  lw=1.5, color='b', ls='--')
plt.xlim([-1.35, 2.7])
plt.xticks([0, 1.35])
plt.errorbar(x=0, y=2020*slope + intercept, yerr=delta[x==2022], lw=1.1, color='grey', capsize=4, fmt = '.', capthick=1.3, ecolor='grey', alpha=1)

plt.ylim([22,29])
# plt.savefig('../../../megafires_data/png/QN_series_with_trend_RR_FET_dec_2022.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()

mu_plus_1sigma = mu + sigma0
mu_plus_2sigma = mu + 2*sigma0

sm = sm.where(sm.time.dt.year != 2020, drop=True)
qn = qn.where(qn.time.dt.year != 2020, drop=True)

fig = plt.figure(figsize=(8,5))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10
plt.scatter(sm, qn, s=12, marker='o', edgecolor='blue', facecolor='none')
plt.scatter(smf.sel(time='2020'), qnf.sel(time='2020'), s=12, marker='s', edgecolor='fuchsia', facecolor='none')
plt.plot(smf, mu, color='red', linewidth = 2)
plt.plot(smf, mu_plus_1sigma, color='red', linewidth = 0.5)
plt.plot(smf, mu_plus_2sigma, color='red', linewidth = 0.5)


plt.ylim([22, 29])
plt.xlim([-0.25, 1.0])
plt.grid(color='grey', lw=0.4, ls='--')
ax = plt.gca()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(direction="in")
plt.xlabel('Global mean surface temperature anomaly (smoothed) [ºC]')
plt.ylabel('December Tmax [ºC]')
#plt.savefig('../../../megafires_data/png/QN_MLE_tasmax_GMST_scatter_fit_1930_2022_dec_2022.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()

tau = 1/norm.sf(qnf.sel(time='2020'), mu.sel(time='2020'), sigma0)
tau2 = 1/norm.sf(qnf.sel(time='2020'), *norm.fit(qnf.values))

smf = gmst.get_gmst_annual_5year_smooth_2022()
mu = mu0 + alpha*smf

tau3 = 1/norm.sf(mu.sel(time='2020'), mu.sel(time='2020'), sigma0)
tau4 = 1/norm.sf(mu.sel(time='2020'), mu.sel(time='1950'), sigma0)

tau5 = 1/norm.sf(qnf.sel(time='2020'), mu.sel(time='2020'), sigma0)
tau6 = 1/norm.sf(qnf.sel(time='2020'), mu.sel(time='1950'), sigma0)

print(tau, tau2, tau3, tau4, tau5, tau6)
