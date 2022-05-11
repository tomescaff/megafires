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
da_jan = da_jan.sel(time=slice('1950', '2021'))
np_jan_all_runs = np.ravel(da_jan.values)

# fit normal dist
normfit = norm.fit(np_jan_all_runs)

# fit gev
gevfit = genextreme.fit(np_jan_all_runs)

# fit gumbel_r
gumrfit = gumbel_r.fit(np_jan_all_runs)

# fit gumbel_l
gumlfit = gumbel_l.fit(np_jan_all_runs)

# test fit
stat, p_norm = kstest(np_jan_all_runs, 'norm', normfit)
stat, p_gev = kstest(np_jan_all_runs, 'genextreme', gevfit)
stat, p_gumr = kstest(np_jan_all_runs, 'gumbel_r', gumrfit)


# compute histogram
hist, bins = np.histogram(np_jan_all_runs, bins=np.arange(20,30.2, 0.2), density=True)

# create figure
fig = plt.figure(figsize=(8,6))

# plot the histogram
width = 0.9 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width, edgecolor='blue', facecolor='lightskyblue', color='blue', alpha = 0.75, label = 'LENS Jan Tmax (all runs)')
plt.xlim([20,30])
# plot the PDF
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)

plt.plot(x, norm.pdf(x, *normfit), 'k', linewidth=2, label = 'Normal fit (p = {:g})'.format(p_norm))
plt.plot(x, genextreme.pdf(x, *gevfit), 'g', linewidth=2, label = 'GEV fit (p = {:g})'.format(p_gev))
plt.plot(x, gumbel_r.pdf(x, *gumrfit), 'r', linewidth=2, label = 'GUM_r fit (p = {:g})'.format(p_gumr))

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set title and labels
plt.xlabel('LENS January Tmax (ºC)')
plt.ylabel('PDF')
plt.savefig('../../../megafires_data/png/LENS_pdf.png', dpi=300)
plt.show()
    

