import sys
import numpy as np
import xarray as xr
import pandas as pd
from scipy.stats import gumbel_r as gum

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.math as pmath

index = ['tau cf', 'tau ac', 'rr c-a', 'far c-a', 'delta c-a']
columns = ['raw', '95ci lower', '95ci upper', '1percentile']
df = pd.DataFrame(columns=columns, index=index)

ac_year = '2017'
cf_year = '1880'

# raw values

smf = gmst.get_gmst_annual_5year_smooth()
qnf = stns.get_QN_tmax_jan()

sm = smf.sel(time=slice('1928','2021'))
qn = qnf.sel(time=slice('1928','2021'))

sm_no2017 = sm.where(sm.time.dt.year != 2017, drop=True)
qn_n02017 = qn.where(qn.time.dt.year != 2017, drop=True)

xopt = pmath.mle_gumbel_2d(qn_n02017.values, sm_no2017.values, [30, 2, 0.5])

mu0, sigma0, alpha = xopt
mu = mu0 + alpha*smf

mu_MLE_ac = mu.sel(time = ac_year)
mu_MLE_cf = mu.sel(time = cf_year)
sigma_MLE = sigma0

# get ev value 
ev = qn.sel(time='2017')

# get return periods
tau_cf = 1/gum.sf(ev, mu_MLE_cf, sigma_MLE)
tau_ac = 1/gum.sf(ev, mu_MLE_ac, sigma_MLE)

# get rr and far
rr_ca = tau_cf/tau_ac
far_ca = (tau_cf-tau_ac)/tau_cf

# get delta
delta = ev - gum.isf(1/tau_ac, mu_MLE_cf, sigma_MLE)

df.loc['tau cf', 'raw'] = tau_cf
df.loc['tau ac', 'raw'] = tau_ac
df.loc['rr c-a', 'raw'] = rr_ca
df.loc['far c-a', 'raw'] = far_ca
df.loc['delta c-a', 'raw'] = delta

# bootstrap MLE
nboot = 1000
filepath = '../../../megafires_data/output/MLE_tasmax_jan_QN_GMST_'+str(nboot)+'.nc'
bspreds = xr.open_dataset(filepath)
bspreds_mu0 = bspreds.mu0.values
bspreds_sigma0 = bspreds.sigma0.values
bspreds_alpha = bspreds.alpha.values

Tac = smf.sel(time = ac_year).values
Tcf = smf.sel(time = cf_year).values

mu_ac_dist = bspreds_mu0 + bspreds_alpha*Tac
mu_cf_dist = bspreds_mu0 + bspreds_alpha*Tcf
sigma_dist = bspreds_sigma0

bspreds_tau_cf = np.zeros((nboot,))
bspreds_tau_ac = np.zeros((nboot,))
bspreds_rr_ca = np.zeros((nboot,))
bspreds_far_ca = np.zeros((nboot,))
bspreds_delta = np.zeros((nboot,))

for i in range(nboot):

    tau_cf_i = 1/gum.sf(ev, mu_cf_dist[i], sigma_dist[i])
    tau_ac_i = 1/gum.sf(ev, mu_ac_dist[i], sigma_dist[i])

    rr_ca_i = tau_cf_i/tau_ac_i
    far_ca_i = (tau_cf_i-tau_ac_i)/tau_cf_i
    delta_i = ev - gum.isf(1/tau_ac_i, mu_cf_dist[i], sigma_dist[i])

    bspreds_tau_cf[i] = tau_cf_i
    bspreds_tau_ac[i] = tau_ac_i
    bspreds_rr_ca[i] = rr_ca_i
    bspreds_far_ca[i] = far_ca_i
    bspreds_delta[i] = delta_i

mapping = [('95ci lower', 0.025), ('95ci upper', 0.975), ('1percentile', 0.01)]

for col, thr in mapping:
    df.loc['tau cf', col] = np.quantile(bspreds_tau_cf, [thr], axis = 0)
    df.loc['tau ac', col] = np.quantile(bspreds_tau_ac, [thr], axis = 0)
    df.loc['rr c-a', col] = np.quantile(bspreds_rr_ca, [thr], axis = 0)
    df.loc['far c-a', col] = np.quantile(bspreds_far_ca, [thr], axis = 0)
    df.loc['delta c-a', col] = np.quantile(bspreds_delta, [thr], axis = 0)

df = df.applymap(lambda x: round(float(x),2))
df.to_csv(f'../../../megafires_data/output/metrics_of_attr_QN_MLE_gumbel_{cf_year}_{ac_year}.csv')
