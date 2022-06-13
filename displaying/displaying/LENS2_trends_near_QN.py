import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats

sys.path.append('../../processing')

import processing.utils as ut
# get LENS trend 1950-2020 jan tmax
ds = xr.open_dataset('../../../megafires_data/LENS2_ALL/LENS2_tmax_mean_mon_1950_2020_jan_trend.nc')
da = ds['slope']*100

lats = [-34.397907, -33.455498, -32.51309]
lons = [-71.25,  -70.0]

# get Quinta Normal time series
qn = ut.get_QN_series()
qn = qn.sel(time=slice('1950','2020'))

y = qn.values
x = qn.time.dt.year.values

slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

fig = plt.figure(figsize=(8,5))

k = 0
for lat in lats:
    for lon in lons:
        real_lon = lon%360
        y = da.sel(lat=lat, lon=real_lon).values
        plt.boxplot(y, notch=False, meanline=True,showmeans=True, positions=[k],  whis=(5,95), patch_artist=True, showfliers=True, labels=[f'lat={lat:.2f}\nlon={lon:.2f}'], boxprops={'facecolor':'grey', 'alpha':0.6}, medianprops=dict(color='blue'), meanprops=dict(color='red'))
        k = k+1
plt.yticks(np.arange(-0.5, 4.5+0.5,0.5))
plt.ylim([-0.7, 4.7])
plt.axhline(slope*100, color='r', lw=1.2, ls='--')
plt.grid(ls='--', color='grey', lw=.4)
plt.savefig('../../../megafires_data/png/LENS2_trends_near_QN.png',dpi=300)
plt.show()
