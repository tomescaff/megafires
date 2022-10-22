import sys
import numpy as np
import xarray as xr
import pandas as pd
from scipy.stats import norm
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.lens as lens
import processing.math as pmath

index = ['tau cf', 'tau ac', 'rr c-a', 'far c-a', 'delta c-a']
columns = ['raw', '95ci lower', '95ci upper', '1percentile']
df = pd.DataFrame(columns=columns, index=index)

ac_year = '2017'

# raw values

lens1_gmst_full = gmst.get_gmst_annual_lens1_ensmean()
lens1_tmax_full = lens.get_LENS_jan_tmax_CU_NN()

lens1_gmst = lens1_gmst_full.sel(time=slice('1920', '2021'))
lens1_tmax = lens1_tmax_full.sel(time=slice('1920', '2021'))

lens1_gmst_arr = np.tile(lens1_gmst.values, lens1_tmax.shape[0])
lens1_tmax_arr = np.ravel(lens1_tmax.values)

xopt = pmath.mle_norm_2d(lens1_tmax_arr, lens1_gmst_arr, [29.55, 1.03, 1.11])

mu0, sigma0, alpha = xopt
mu = mu0 + alpha*lens1_gmst_full

mu_MLE_ac = mu.sel(time = ac_year)
sigma_MLE = sigma0

# counterfactual
cr = lens.get_LENS_jan_tmax_CU_control_run()
cr_normfit = norm.fit(cr.values)

# define tau
tau_ac = 241

# get ev value 
ev = norm.isf(1/tau_ac, mu_MLE_ac, sigma_MLE)

# get return periods
tau_cf = 1/norm.sf(ev, *cr_normfit)
tau_ac = 1/norm.sf(ev, mu_MLE_ac, sigma_MLE)

# get rr and far
rr_ca = tau_cf/tau_ac
far_ca = (tau_cf-tau_ac)/tau_cf

# get delta
delta = ev - norm.isf(1/tau_ac, *cr_normfit)

df.loc['tau cf', 'raw'] = tau_cf
df.loc['tau ac', 'raw'] = tau_ac
df.loc['rr c-a', 'raw'] = rr_ca
df.loc['far c-a', 'raw'] = far_ca
df.loc['delta c-a', 'raw'] = delta

# bootstrap MLE
nboot = 1000
filepath = '../../../megafires_data/output/MLE_tasmax_jan_LENS1_GMST_'+str(nboot)+'_normal_evaluation_CU_NN.nc'
bspreds = xr.open_dataset(filepath)
bspreds_mu0 = bspreds.mu0.values
bspreds_sigma0 = bspreds.sigma0.values
bspreds_alpha = bspreds.alpha.values

Tac = lens1_gmst_full.sel(time = ac_year).values

mu_ac_dist = bspreds_mu0 + bspreds_alpha*Tac
sigma_dist = bspreds_sigma0

bspreds_tau_cf = np.zeros((nboot,))
bspreds_tau_ac = np.zeros((nboot,))
bspreds_rr_ca = np.zeros((nboot,))
bspreds_far_ca = np.zeros((nboot,))
bspreds_delta = np.zeros((nboot,))

for i in range(nboot):

    cr_i = bootstrap(cr.values)
    cr_normfit_i = norm.fit(cr_i)
    tau_cf_i = 1/norm.sf(ev, *cr_normfit_i)
    tau_ac_i = 1/norm.sf(ev, mu_ac_dist[i], sigma_dist[i])

    rr_ca_i = tau_cf_i/tau_ac_i
    far_ca_i = (tau_cf_i-tau_ac_i)/tau_cf_i

    # ev_i = norm.isf(1/tau_ac, mu_ac_dist[i], sigma_dist[i])
    delta_i = ev - norm.isf(1/tau_ac, *cr_normfit_i)

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
df.to_csv(f'../../../megafires_data/output/metrics_LENS1_MLE_crun_{ac_year}_normfit_1920_2021_by_return_period_CU_NN.csv')