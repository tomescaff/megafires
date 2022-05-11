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

# get Quinta Normal time series
basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan = da[:, ::12]
da_jan = da_jan.sel(time=slice('1920', '1950'))
np_jan_all_runs = np.ravel(da_jan.values)

# get Quinta Normal time series
basedir = '../../../megafires_data/LENS/'
filename = 'LENS_tmax_mon_QN_control_run.nc'
filepath = basedir + filename
ds = xr.open_dataset(filepath)
da = ds['TREFMXAV']-273.15
da_jan = da[::12]
np_jan_control_run = np.ravel(da_jan.values)

# fit normal dist
normfit_ar = norm.fit(np_jan_all_runs)
normfit_cr = norm.fit(np_jan_control_run)

# create figure
fig = plt.figure(figsize=(8,6))

plt.xlim([20,30])
# plot the PDF
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)

plt.plot(x, norm.pdf(x, *normfit_ar), 'b', linewidth=2, label = f'Historical (mu={normfit_ar[0]:.2f}, sig={normfit_ar[1]:.2f})')
plt.plot(x, norm.pdf(x, *normfit_cr), 'r', linewidth=2, label = f'Counterfactual (mu={normfit_cr[0]:.2f}, sig={normfit_cr[1]:.2f})')

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=1)

# set title and labels
plt.xlabel('LENS January Tmax (ºC)')
plt.ylabel('PDF')
plt.savefig('../../../megafires_data/png/LENS_pdf_all_runs_vs_control_run_1920_1950.png', dpi=300)
plt.show()
    

