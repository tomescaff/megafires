import sys
import numpy as np
import xarray as xr
import pandas as pd
from scipy.stats import norm
from scipy.stats import genextreme as gev

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.math as pmath

index = ['tau cf', 'tau ac', 'rr c-a', 'far c-a', 'delta c-a']
columns = ['raw', '95ci lower', '95ci upper', '1percentile']
df = pd.DataFrame(columns=columns, index=index)

ac_year = '2017'
cf_year = '1943'

# raw values

smf = gmst.get_gmst_annual_5year_smooth()
qnf = stns.get_QN_tmax_jan()

sm = smf.sel(time=slice('1928','2021'))
qn = qnf.sel(time=slice('1928','2021'))

sm_no2017 = sm.where(sm.time.dt.year != 2017, drop=True)
qn_n02017 = qn.where(qn.time.dt.year != 2017, drop=True)

xopt = pmath.mle_gev_2d(qn_n02017.values, sm_no2017.values, [15, 2, 0.5, 0])

mu0, sigma0, alpha, eta = xopt
mu = mu0 + alpha*smf

mu_MLE_ac = mu.sel(time = ac_year)
mu_MLE_cf = mu.sel(time = cf_year)
sigma_MLE = sigma0

# get ev value 
ev = qn.sel(time='2017')

# get return periods
tau_cf = 1/gev.sf(ev, eta, mu_MLE_cf, sigma_MLE)
tau_ac = 1/gev.sf(ev, eta, mu_MLE_ac, sigma_MLE)

# get rr and far
rr_ca = tau_cf/tau_ac
far_ca = (tau_cf-tau_ac)/tau_cf

# get delta
delta = ev - gev.isf(1/tau_ac, eta, mu_MLE_cf, sigma_MLE)

df.loc['tau cf', 'raw'] = tau_cf
df.loc['tau ac', 'raw'] = tau_ac
df.loc['rr c-a', 'raw'] = rr_ca
df.loc['far c-a', 'raw'] = far_ca
df.loc['delta c-a', 'raw'] = delta

# after 33.27638225 pdf is 0 -> taus are inf
