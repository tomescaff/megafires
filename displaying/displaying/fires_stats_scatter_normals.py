import sys
from syslog import LOG_SYSLOG
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
from scipy.stats import pearsonr
sys.path.append('../../processing')

import processing.utils as ut
import processing.stations as stns
import processing.fires as fires
from scipy.stats import multivariate_normal
from scipy.stats import linregress
from scipy.stats import norm

# possible x
qn_jan = stns.get_QN_tmax_jan()
cu_jan = stns.get_CU_tmax_jan()
ch_jan = stns.get_CH_tmax_jan()
qn_feb = stns.get_QN_tmax_feb()
cu_feb = stns.get_CU_tmax_feb()
ch_feb = stns.get_CH_tmax_feb()
regind_jan = stns.get_regindex_tmax_jan()
regind_feb = stns.get_regindex_tmax_feb()
qn_janfeb = stns.get_QN_tmax_janfeb()

# possible y
ff_seas = fires.get_fires_frequency_by_season()
ba_seas = fires.get_burned_area_by_season()
ba_jan = fires.get_burned_area_january()
ba_feb = fires.get_burned_area_february()

# select x and y
x_full = qn_jan
y_full = ba_seas

# remove 2017
remove_2017 = True
if remove_2017:
    x = x_full.where(x_full.time.dt.year != 2017, drop=True)
    y = y_full.where(y_full.time.dt.year != 2017, drop=True)
else:
    x = x_full
    y = y_full

x_2017 = x_full.where(x_full.time.dt.year == 2017, drop=True)
y_2017 = y_full.where(y_full.time.dt.year == 2017, drop=True)

# slice time
init = '1964'
end = '2021'

x = x.sel(time=slice(init, end))
y = y.sel(time=slice(init, end))

print(pearsonr(x,y))
logy=True
if logy:
    y = np.log10(y)
    y_2017 = np.log10(y_2017)
print(pearsonr(x,y))

x_min = 27.2 # x.min().values - x.std().values
x_max = 34 # x.max().values + x.std().values
y_min = 3.9-1.5# y.min().values - y.std().values
y_max = 5.9# y.max().values + y.std().values

# display contour

# x_, y_ = np.mgrid[x_min:x_max:.01, y_min:y_max:.01]
# x_, y_ = np.meshgrid(np.linspace(x_min,x_max, 100), np.linspace(y_min, y_max, 100))
# pos = np.dstack((x_, y_))
# mu_x = np.mean(x)
# mu_y = np.mean(y)
# cov = np.cov(x, y)
# rv = multivariate_normal([mu_x, mu_y], cov)

X = np.vstack((x.values,y.values)).T
sigma = np.cov(X.T)
det = np.linalg.det(sigma)
inv = np.linalg.inv(sigma)
def fun(x0,x1):
    x = np.array([x0,x1])
    const=1/(2*np.pi*det)**0.5
    arg = -(x-X.mean(axis=0)).T@inv@(x-X.mean(axis=0))
    return np.exp(arg)*const

def d2(x0,x1):
    x = np.array([x0,x1])
    return (x-X.mean(axis=0)).T@inv@(x-X.mean(axis=0))
funvec = np.vectorize(fun)
d2vec = np.vectorize(d2)
x_ = np.linspace(x_min,x_max,100)
y_ = np.linspace(y_min,y_max,100)

xx,yy = np.meshgrid(x_,y_)
zz = funvec(xx,yy)
dd = d2vec(xx,yy)
    
# display regress

slope, intercept, rvalue, pvalue, stderr  = linregress(x,y)
p = np.polyfit(x,y,1)
xval = np.linspace(x_min, x_max, 100)
yval = slope*xval + intercept
yval = p[0]*xval + p[1]

# figure
fig = plt.figure(figsize=(10,7.5))

# Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
# the size of the marginal axes and the main axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(2, 2,  width_ratios=(7.5, 1.5), height_ratios=(2, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.07, hspace=0.07)

ax_scatt = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_scatt)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_scatt)

plt.sca(ax_scatt)
plt.contourf(xx,yy,dd,cmap='Greys',levels=[0,4.61],alpha=0.2)
plt.contour(xx,yy,dd,colors=['grey'],levels=[9.21], linewidths=[0.5], linestyles='--')
plt.scatter(x,y, s=50, edgecolor='blue', facecolor='lightskyblue')
plt.axvline(x.mean(), color='b', lw=0.5, ls='--')
plt.axhline(y.mean(), color='b', lw=0.5, ls='--')
plt.plot(xval, yval, color='b', lw=1.5, ls='--')
plt.axvline(x_2017, color='r', lw=0.5, ls='--')
plt.axhline(y_2017, color='r', lw=0.5, ls='--')
plt.scatter(x_2017,y_2017, s=50, edgecolor='red', facecolor='coral', zorder=4)
plt.xticks(np.arange(27.2, 35, 1))
plt.yticks(np.arange(3.50, 6.25, 0.25), ["{:.2E}".format(10**x) for x in np.arange(3.50, 6.25, 0.25)])
plt.xlim([27.2, 34])
plt.ylim([3.6, 6.0])
# Hide the right and top spines
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
plt.xlabel('Tmax January (ºC)')
plt.ylabel('Burned area per season (Ha)')

for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(9) 
        tick.label.set_weight('light') 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(9) 
    tick.label.set_weight('light')

# top
plt.sca(ax_histx)
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(0) 
    tick.label.set_color('white') 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(0) 
    tick.label.set_color('white') 
# fit normal dist
normfit = norm.fit(x)
# plot the PDF
xval = np.linspace(x_min, x_max, 100)
plt.fill_between(xval, xval*0, norm.pdf(xval, *normfit), facecolor='grey', linewidth=2, alpha=0.2)
plt.ylim([0, norm.pdf(xval, *normfit).max()])
plt.ylabel('PDF')
plt.axvline(x_2017, color='r', lw=0.5, ls='--')
plt.axvline(x.mean(), color='b', lw=0.5, ls='--')

# right
plt.sca(ax_histy)
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(0)
    tick.label.set_color('white') 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(0)
    tick.label.set_color('white') 
# fit normal dist
normfit = norm.fit(y)
# plot the PDF
xval = np.linspace(y_min-1, y_max, 100)
plt.fill_betweenx(  xval, xval*0, norm.pdf(xval, *normfit), facecolor='grey', linewidth=2, alpha=0.2)
plt.xlim([0, norm.pdf(xval, *normfit).max()])
plt.xlabel('PDF')
plt.axhline(y_2017, color='r', lw=0.5, ls='--')
plt.axhline(y.mean(), color='b', lw=0.5, ls='--')
plt.savefig('../../../megafires_data/png/fires_stats_scatter_normals.png', dpi=300)
plt.show()

