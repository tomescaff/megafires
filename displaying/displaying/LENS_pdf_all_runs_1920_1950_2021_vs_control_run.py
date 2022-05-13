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
da_jan = da[:, ::12]

# get 1920-1950 period
da_jan_1920_1950 = da_jan.sel(time=slice('1920', '1950'))
np_jan_all_runs_1920_1950 = np.ravel(da_jan_1920_1950.values)

# get 1950-2021 period
da_jan_1950_2021 = da_jan.sel(time=slice('1950', '2021'))
np_jan_all_runs_1950_2021 = np.ravel(da_jan_1950_2021.values)

# get control run time series
basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN_control_run.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan = da[::12]
np_jan_control_run = np.ravel(da_jan.values)

# fit normal dist
normfit_ar_1920_1950 = norm.fit(np_jan_all_runs_1920_1950)
normfit_ar_1950_2021 = norm.fit(np_jan_all_runs_1950_2021)
normfit_cr = norm.fit(np_jan_control_run)

# create figure
fig = plt.figure(figsize=(8,6))

plt.xlim([20,30])
# plot the PDF
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)

plt.plot(x, norm.pdf(x, *normfit_ar_1950_2021), 'b', linewidth=2, label = 'LENS (40 runs) 1950-2021')
plt.plot(x, norm.pdf(x, *normfit_ar_1920_1950), 'k', linewidth=2, label = 'LENS (40 runs) 1920-1950')
plt.plot(x, norm.pdf(x, *normfit_cr), 'r', linewidth=2, label = 'Control run')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set title and labels
plt.xlabel('LENS January Tmax (ºC)')
plt.ylabel('PDF')
plt.savefig('../../../megafires_data/png/LENS_pdf_all_runs_1920_1950_2021_vs_control_run.png', dpi=300)
plt.show()
    

