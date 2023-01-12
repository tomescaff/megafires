import sys
import numpy as np
import xarray as xr
import pandas as pd
from scipy.stats import norm

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.lens as lens
import processing.math as pmath

index = ['tau cf', 'tau ac', 'rr c-a', 'far c-a', 'delta c-a']
columns = ['raw', '95ci lower', '95ci upper', '1percentile']
df = pd.DataFrame(columns=columns, index=index)

ac_year = '2017'
cf_year = '1880'

# raw values

lens2_gmst_full = gmst.get_gmst_annual_lens2_ensmean()
lens2_tmax_full = lens.get_LENS2_jan_tmax_QNW()

lens2_gmst = lens2_gmst_full.sel(time=slice('1850', '2021'))
lens2_tmax = lens2_tmax_full.sel(time=slice('1850', '2021'))

lens2_gmst_arr = np.tile(lens2_gmst.values, lens2_tmax.shape[0])
lens2_tmax_arr = np.ravel(lens2_tmax.values)

xopt = pmath.mle_norm_2d(lens2_tmax_arr, lens2_gmst_arr, [29.55, 1.03, 1.11])

mu0, sigma0, alpha = xopt
mu = mu0 + alpha*lens2_gmst_full

mu_MLE_ac = mu.sel(time = ac_year)
mu_MLE_cf = mu.sel(time = cf_year)
sigma_MLE_acf = sigma0

# define tau
tau_ac = 1397

# get ev value 
ev = norm.isf(1/tau_ac, mu_MLE_ac, sigma_MLE_acf)

#################

ac_year = '2017'
fu_year = '2070'

# raw values

lens2_gmst_full = gmst.get_gmst_annual_lens2_ensmean()
lens2_tmax_full = lens.get_LENS2_jan_tmax_QNW()

lens2_gmst = lens2_gmst_full.sel(time=slice('2010', '2100'))
lens2_tmax = lens2_tmax_full.sel(time=slice('2010', '2100'))

lens2_gmst_arr = np.tile(lens2_gmst.values, lens2_tmax.shape[0])
lens2_tmax_arr = np.ravel(lens2_tmax.values)

xopt = pmath.mle_norm_2d(lens2_tmax_arr, lens2_gmst_arr, [29.55, 1.03, 1.11])

mu0, sigma0, alpha = xopt
mu = mu0 + alpha*lens2_gmst_full

mu_MLE_fu = mu.sel(time = fu_year)
sigma_MLE_fu = sigma0

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
# create figure
fig = plt.figure(figsize=(9,6))

plt.xlim([22,36])
# plot the PDF
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)

plt.plot(x, norm.pdf(x, mu_MLE_cf, sigma_MLE_acf), 'b', linewidth=2)
plt.plot(x, norm.pdf(x, mu_MLE_ac, sigma_MLE_acf), 'k', linewidth=2)
plt.plot(x, norm.pdf(x, mu_MLE_fu, sigma_MLE_fu), 'r', linewidth=2)

plt.axvline(ev, color='fuchsia', linewidth=1.5)

plt.ylim([0,0.33])
# set grid
plt.grid(lw=0.2, ls='--', color='grey')
ax = plt.gca()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(direction="in")
# set title and labels
plt.savefig('../../../megafires_data/png/doca_pdf_present_future_normal.png', dpi=300)
plt.show()