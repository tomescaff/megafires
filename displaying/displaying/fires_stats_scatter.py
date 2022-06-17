import sys
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
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
x = qn_jan
y = ba_seas

# remove 2017
remove_2017 = True
if remove_2017:
    x = x.where(x.time.dt.year != 2017, drop=True)
    y = y.where(y.time.dt.year != 2017, drop=True)

# slice time
init = '1964'
end = '2021'

x = x.sel(time=slice(init, end))
y = y.sel(time=slice(init, end))

print(pearsonr(x,y))
log=True
if log:
    y = np.log10(y)
print(pearsonr(x,y))

x_min = 28 # x.min().values - x.std().values
x_max = 32 # x.max().values + x.std().values
y_min = 3.9# y.min().values - y.std().values
y_max = 5.2# y.max().values + y.std().values

# display contour

x_, y_ = np.mgrid[x_min:x_max:.01, y_min:y_max:.01]
x_, y_ = np.meshgrid(np.linspace(x_min,x_max, 100), np.linspace(y_min, y_max, 100))
pos = np.dstack((x_, y_))
mu_x = np.mean(x)
mu_y = np.mean(y)
cov = np.cov(x, y)
rv = multivariate_normal([mu_x, mu_y], cov)

# display regress

slope, intercept, rvalue, pvalue, stderr  = linregress(x,y)
p = np.polyfit(x,y,1)
xval = np.linspace(x_min, x_max, 100)
yval = slope*xval + intercept
yval = p[0]*xval + p[1]

# figure
fig = plt.figure(figsize=(12,8))
plt.contour(x_, y_, rv.pdf(pos), colors='k', linewidths=0.5, linestyles='--')
plt.scatter(x,y, s=50, edgecolor='blue', facecolor='lightskyblue')
plt.axvline(mu_x, color='r', lw=0.5, ls='--')
plt.axhline(mu_y, color='r', lw=0.5, ls='--')
plt.plot(xval, yval, color='b', lw=0.5, ls='--')
plt.show()
