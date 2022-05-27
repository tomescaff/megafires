import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
from scipy.stats import norm
from scipy.stats import genextreme
from scipy.stats import gumbel_r, gumbel_l
from scipy import stats
from scipy.stats import kstest
import numpy as np
import xarray as xr

# get LENS time series
basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
lens = da[:, ::12]

# fixing LENS series
qn = ut.get_QN_series()
lens_mean_1950_1970 = lens.sel(time=slice('1950', '1970')).mean('time')
qn_mean_1950_1970 = qn.sel(time=slice('1950', '1970')).mean('time')
lens_anom = lens - lens_mean_1950_1970
da_jan = lens_anom + qn_mean_1950_1970

# get LENS 1920-1950 period
da_jan_1920_1950 = da_jan.sel(time=slice('1920', '1950'))
np_jan_all_runs_1920_1950 = np.ravel(da_jan_1920_1950.values)

# get LENS 1990-2020 period
da_jan_1990_2020 = da_jan.sel(time=slice('1990', '2020'))
np_jan_all_runs_1990_2020 = np.ravel(da_jan_1990_2020.values)

# get QN 1990-2020 period
qn_1990_2020 = qn.sel(time=slice('1990', '2020'))
np_jan_qn_1990_2020 = np.ravel(qn_1990_2020.values)

# get LENS 2070-2100 period
da_jan_2070_2100 = da_jan.sel(time=slice('2070', '2100'))
np_jan_all_runs_2070_2100 = np.ravel(da_jan_2070_2100.values)

# fit normal
normfit_ar_1920_1950 = norm.fit(np_jan_all_runs_1920_1950)
normfit_ar_1990_2020 = norm.fit(np_jan_all_runs_1990_2020)
normfit_qn_1990_2020 = norm.fit(np_jan_qn_1990_2020)
normfit_ar_2070_2100 = norm.fit(np_jan_all_runs_2070_2100)


# create figure
fig = plt.figure(figsize=(9,6))

plt.xlim([25,38])
# plot the PDF
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)

plt.plot(x, norm.pdf(x, *normfit_ar_1920_1950), 'b', linewidth=2, label = 'LENS (40 runs) 1920-1950')
plt.plot(x, norm.pdf(x, *normfit_ar_1990_2020), 'k', linewidth=2, label = 'LENS (40 runs) 1990-2020')
plt.plot(x, norm.pdf(x, *normfit_qn_1990_2020), 'green', linewidth=2, label = 'QN 1990-2020')
plt.plot(x, norm.pdf(x, *normfit_ar_2070_2100), 'r', linewidth=2, label = 'LENS (40 runs) 2070-2100')

plt.axvline(qn.sel(time='2017'), color='grey', linewidth=1.2, label = 'Temp 2017')

plt.ylim([0,0.5])
# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=3)

# set title and labels
plt.xlabel('LENS January Tmax (ºC)')
plt.ylabel('PDF')
plt.savefig('../../../megafires_data/png/LENS_corrected_pdf_past_present_future_QN_normal.png', dpi=300)
plt.show()
    

