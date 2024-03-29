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
remove_2017 = False
if remove_2017:
    x = x_full.where(x_full.time.dt.year != 2017, drop=True)
    y = y_full.where(y_full.time.dt.year != 2017, drop=True)
else:
    x = x_full
    y = y_full

x_2017 = x_full.where(x_full.time.dt.year == 2017, drop=True)
y_2017 = y_full.where(y_full.time.dt.year == 2017, drop=True)

# slice time
init = '1995'
end = '2021'

x = x.sel(time=slice(init, end))
y = y.sel(time=slice(init, end))

print(pearsonr(x,y))
logy=True
if logy:
    y = np.log10(y)
    y_2017 = np.log10(y_2017)
print(pearsonr(x,y))

x_min = 27.5 # x.min().values - x.std().values
x_max = 34 # x.max().values + x.std().values
y_min = 3.9# y.min().values - y.std().values
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
fig = plt.figure(figsize=(10,7))
#plt.contourf(xx,yy,zz,15,cmap='Greys',alpha=0.8)
plt.contourf(xx,yy,dd,cmap='Greys',levels=[0,4.61],alpha=0.2)
#plt.contour(x_, y_, rv.pdf(pos), colors='k', linewidths=0.5, linestyles='--')
plt.scatter(x,y, s=50, edgecolor='blue', facecolor='lightskyblue')
plt.axvline(x.mean(), color='b', lw=0.5, ls='--')
plt.axhline(y.mean(), color='b', lw=0.5, ls='--')
plt.plot(xval, yval, color='b', lw=1.5, ls='--')
plt.axvline(x_2017, color='r', lw=0.5, ls='--')
plt.axhline(y_2017, color='r', lw=0.5, ls='--')
plt.scatter(x_2017,y_2017, s=50, edgecolor='red', facecolor='coral', zorder=4)
plt.xticks(np.arange(28, 35, 1))
plt.yticks(np.arange(4.00, 6.25, 0.25), ["{:.2E}".format(10**x) for x in np.arange(4.00, 6.25, 0.25)])
plt.xlim([27.5, 34])
plt.ylim([3.8, 6.0])
# Hide the right and top spines
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
plt.xlabel('Tmax Jan (ºC)')
plt.ylabel('Burned area per season (Ha)')

for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(9) 
        tick.label.set_weight('light') 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(9) 
    tick.label.set_weight('light')
str = '2017_used' if not remove_2017 else '2017_removed'
plt.savefig(f'../../../megafires_data/png/fires_stats_scatter_{init}_{end}_{str}.png')
plt.show()
