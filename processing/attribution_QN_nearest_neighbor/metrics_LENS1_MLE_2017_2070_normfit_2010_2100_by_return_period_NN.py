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

index = ['tau ac', 'tau fu', 'rr a-f', 'far a-f', 'delta a-f']
columns = ['raw', '95ci lower', '95ci upper', '1percentile']
df = pd.DataFrame(columns=columns, index=index)

ac_year = '2017'
fu_year = '2070'

# raw values

lens1_gmst_full = gmst.get_gmst_annual_lens1_ensmean()
lens1_tmax_full = lens.get_LENS_jan_tmax_QNW()

lens1_gmst = lens1_gmst_full.sel(time=slice('2010', '2100'))
lens1_tmax = lens1_tmax_full.sel(time=slice('2010', '2100'))

lens1_gmst_arr = np.tile(lens1_gmst.values, lens1_tmax.shape[0])
lens1_tmax_arr = np.ravel(lens1_tmax.values)

xopt = pmath.mle_norm_2d(lens1_tmax_arr, lens1_gmst_arr, [29.55, 1.03, 1.11])

mu0, sigma0, alpha = xopt
mu = mu0 + alpha*lens1_gmst_full

mu_MLE_ac = mu.sel(time = ac_year)
mu_MLE_fu = mu.sel(time = fu_year)
sigma_MLE = sigma0

# define tau
tau_ac = 1397

# get ev value 
ev = norm.isf(1/tau_ac, mu_MLE_ac, sigma_MLE)

# get return periods
tau_fu = 1/norm.sf(ev, mu_MLE_fu, sigma_MLE)
tau_ac = 1/norm.sf(ev, mu_MLE_ac, sigma_MLE)

# get rr and far
rr_af = tau_ac/tau_fu
far_af = (tau_ac-tau_fu)/tau_ac

# get delta
delta = ev - norm.isf(1/tau_ac, mu_MLE_fu, sigma_MLE)

df.loc['tau fu', 'raw'] = tau_fu
df.loc['tau ac', 'raw'] = tau_ac
df.loc['rr a-f', 'raw'] = rr_af
df.loc['far a-f', 'raw'] = far_af
df.loc['delta a-f', 'raw'] = delta

# bootstrap MLE
nboot = 1000
filepath = '../../../megafires_data/output/MLE_tasmax_jan_LENS1_GMST_'+str(nboot)+'_normal_future_QN_NN.nc'
bspreds = xr.open_dataset(filepath)
bspreds_mu0 = bspreds.mu0.values
bspreds_sigma0 = bspreds.sigma0.values
bspreds_alpha = bspreds.alpha.values

Tac = lens1_gmst_full.sel(time = ac_year).values
Tfu = lens1_gmst_full.sel(time = fu_year).values

mu_ac_dist = bspreds_mu0 + bspreds_alpha*Tac
mu_fu_dist = bspreds_mu0 + bspreds_alpha*Tfu
sigma_dist = bspreds_sigma0

bspreds_tau_fu = np.zeros((nboot,))
bspreds_tau_ac = np.zeros((nboot,))
bspreds_rr_af = np.zeros((nboot,))
bspreds_far_af = np.zeros((nboot,))
bspreds_delta = np.zeros((nboot,))

for i in range(nboot):

    tau_fu_i = 1/norm.sf(ev, mu_fu_dist[i], sigma_dist[i])
    tau_ac_i = 1/norm.sf(ev, mu_ac_dist[i], sigma_dist[i])

    rr_af_i = tau_ac_i/tau_fu_i
    far_af_i = (tau_ac_i-tau_fu_i)/tau_ac_i
    delta_i = ev - norm.isf(1/tau_ac_i, mu_fu_dist[i], sigma_dist[i])

    bspreds_tau_fu[i] = tau_fu_i
    bspreds_tau_ac[i] = tau_ac_i
    bspreds_rr_af[i] = rr_af_i
    bspreds_far_af[i] = far_af_i
    bspreds_delta[i] = delta_i

mapping = [('95ci lower', 0.025), ('95ci upper', 0.975), ('1percentile', 0.01)]

for col, thr in mapping:
    df.loc['tau fu', col] = np.quantile(bspreds_tau_fu, [thr], axis = 0)
    df.loc['tau ac', col] = np.quantile(bspreds_tau_ac, [thr], axis = 0)
    df.loc['rr a-f', col] = np.quantile(bspreds_rr_af, [thr], axis = 0)
    df.loc['far a-f', col] = np.quantile(bspreds_far_af, [thr], axis = 0)
    df.loc['delta a-f', col] = np.quantile(bspreds_delta, [thr], axis = 0)

df = df.applymap(lambda x: round(float(x),2))
df.to_csv(f'../../../megafires_data/output/metrics_LENS1_MLE_{ac_year}_{fu_year}_normfit_2010_2100_by_return_period_QN_NN.csv')