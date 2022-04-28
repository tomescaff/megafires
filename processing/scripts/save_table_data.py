import sys
import pandas as pd
from xarray import DataArray
import xarray as xr

sys.path.append('..')

import processing.utils as ut
import numpy as np
from scipy.stats import norm
from scipy.stats import genextreme
from scipy.stats import gumbel_r

df = pd.DataFrame( columns=['Full series', '2017 removed'])

# get series
qn = ut.get_QN_series()
qn_2017r = qn.where(qn.time.dt.year != 2017, drop = True)

# compute clims
qn_clim_1950_2021 = qn.mean('time')
qn_clim_1981_2010 = qn.sel(time=slice('1981','2010')).mean('time')
df.loc['QN Clim 1950-2021','Full series'] = qn_clim_1950_2021.values.round(2)
df.loc['QN Clim 1981-2010','Full series'] = qn_clim_1981_2010.values.round(2)


# compute anoms
qn_anom_1950_2021 = qn.sel(time='2017') - qn_clim_1950_2021
qn_anom_1981_2010 = qn.sel(time='2017') - qn_clim_1981_2010
df.loc['QN 2017 anom ref. 1950-2021','Full series'] = qn_anom_1950_2021.squeeze().values.round(2)
df.loc['QN 2017 anom ref. 1981-2010','Full series'] = qn_anom_1981_2010.squeeze().values.round(2)

# compute taus
u, tau = ut.get_QN_tau()
qn_tau_non_parametric = tau[-1]
df.loc['QN 2017 non parametric tau','Full series'] = qn_tau_non_parametric.round(2)

# fit normal dist
normfit = norm.fit(qn.values)
qn_tau_parametric_norm = 1/norm.sf(qn.sel(time='2017'), *normfit)

# fit gev
gevfit = genextreme.fit(qn.values)
qn_tau_parametric_gev = 1/genextreme.sf(qn.sel(time='2017'), *gevfit)

# fit gumbel_r
gumrfit = gumbel_r.fit(qn.values)
qn_tau_parametric_gum = 1/gumbel_r.sf(qn.sel(time='2017'), *gumrfit)

df.loc['QN 2017 parametric tau (Norm)','Full series'] = qn_tau_parametric_norm.squeeze().round(2)
df.loc['QN 2017 parametric tau (GEV)','Full series'] = qn_tau_parametric_gev.squeeze().round(2)
df.loc['QN 2017 parametric tau (Gumbel)','Full series'] = qn_tau_parametric_gum.squeeze().round(2)

#################################
# compute taus with 2017 removed
#################################

# fit normal dist 2017 removed
normfit = norm.fit(qn_2017r.values)
qn_tau_parametric_norm_2017r = 1/norm.sf(qn.sel(time='2017'), *normfit)

# fit gev 2017 removed
gevfit = genextreme.fit(qn_2017r.values)
qn_tau_parametric_gev_2017r = 1/genextreme.sf(qn.sel(time='2017'), *gevfit)

# fit gumbel_r 2017 removed
gumrfit = gumbel_r.fit(qn_2017r.values)
qn_tau_parametric_gum_2017r = 1/gumbel_r.sf(qn.sel(time='2017'), *gumrfit)

df.loc['QN 2017 parametric tau (Norm)','2017 removed'] = qn_tau_parametric_norm_2017r.squeeze().round(2)
df.loc['QN 2017 parametric tau (GEV)','2017 removed'] = qn_tau_parametric_gev_2017r.squeeze().round(2)
df.loc['QN 2017 parametric tau (Gumbel)','2017 removed'] = qn_tau_parametric_gum_2017r.squeeze().round(2)

# get linear trend parameters from full series
dtrend = ut.get_linear_trend()

# get predicted values from full series
qn_pred = xr.DataArray(dtrend['y_pred'], coords=[qn.time], dims=['time'])

# get detrended series
qn_detrended = qn - qn_pred + qn.mean('time')

# compute clims detrended
df.loc['QN Clim 1950-2021','Full series linear detrended'] = qn_detrended.mean('time').values.round(2)
df.loc['QN Clim 1981-2010','Full series linear detrended'] = qn_detrended.sel(time=slice('1981','2010')).mean('time').values.round(2)

# compute anoms detrended
qn_detrended_anom_1950_2021 = qn_detrended.sel(time='2017') - qn_detrended.mean('time')
qn_detrended_anom_1981_2010 = qn_detrended.sel(time='2017') - qn_detrended.sel(time=slice('1981','2010')).mean('time')
df.loc['QN 2017 anom ref. 1950-2021','Full series linear detrended'] = qn_detrended_anom_1950_2021.squeeze().values.round(2)
df.loc['QN 2017 anom ref. 1981-2010','Full series linear detrended'] = qn_detrended_anom_1981_2010.squeeze().values.round(2)

# get 2017 values
df.loc['QN 2017 value', 'Full series'] = qn.sel(time='2017').squeeze().values.round(2)
df.loc['QN 2017 value', 'Full series linear detrended'] = qn_detrended.sel(time='2017').squeeze().values.round(2)

# get 2017 diff: full series - x
df.loc['QN 2017 Full series diff', 'Full series'] = 0.0
df.loc['QN 2017 Full series diff', 'Full series linear detrended'] = (df.loc['QN 2017 value', 'Full series'] - df.loc['QN 2017 value', 'Full series linear detrended']).round(2)

# get max values 
df.loc['QN max value', 'Full series'] = qn.max().values.round(2)
df.loc['QN max value', '2017 removed'] = qn_2017r.max().values.round(2)
df.loc['QN max value', 'Full series linear detrended'] = qn_detrended.max().values.round(2)

# get argmax values
df.loc['QN argmax value', 'Full series'] = qn.time.dt.year[qn.argmax()].values
df.loc['QN argmax value', '2017 removed'] = qn_2017r.time.dt.year[qn_2017r.argmax()].values
df.loc['QN argmax value', 'Full series linear detrended'] = qn_detrended.time.dt.year[qn_detrended.argmax()].values

# get trend parameters
df.loc['QN 2017 predict', 'Full series linear model'] = qn_pred.sel(time='2017').squeeze().values.round(2)
df.loc['slope (ºC/year)', 'Full series linear model'] = dtrend['b'].round(2)
df.loc['intercept (year)', 'Full series linear model'] = dtrend['a'].round(2)
df.loc['R2', 'Full series linear model'] = (dtrend['r']**2).round(2)

# compute taus detrended
u_det, tau_det = ut.get_QN_tau_detrended()
df.loc['QN 2017 non parametric tau','Full series linear detrended'] = tau_det[-1].round(2)

# fit normal dist
normfit = norm.fit(qn_detrended.values)
qn_detrended_tau_parametric_norm = 1/norm.sf(qn.sel(time='2017'), *normfit)

# fit gev
gevfit = genextreme.fit(qn_detrended.values)
qn_detrended_tau_parametric_gev = 1/genextreme.sf(qn.sel(time='2017'), *gevfit)

# fit gumbel_r
gumrfit = gumbel_r.fit(qn_detrended.values)
qn_detrended_tau_parametric_gum = 1/gumbel_r.sf(qn.sel(time='2017'), *gumrfit)

df.loc['QN 2017 parametric tau (Norm)','Full series linear detrended'] = qn_detrended_tau_parametric_norm.squeeze().round(2)
df.loc['QN 2017 parametric tau (GEV)','Full series linear detrended'] = qn_detrended_tau_parametric_gev.squeeze().round(2)
df.loc['QN 2017 parametric tau (Gumbel)','Full series linear detrended'] = qn_detrended_tau_parametric_gum.squeeze().round(2)

############################################################
# get linear trend parameters from series with 2017 removed
############################################################
dtrend_2017r = ut.get_linear_trend_2017r()

# get predicted values from series with 2017 removed
qn_pred_2017r = xr.DataArray(dtrend_2017r['y_pred'], coords=[qn_2017r.time], dims=['time'])

# get detrended series from series with 2017 removed
qn_detrended_2017r = qn_2017r - qn_pred_2017r + qn_2017r.mean('time')

# fit normal dist
normfit = norm.fit(qn_detrended_2017r.values)
qn_detrended_2017r_tau_parametric_norm = 1/norm.sf(qn.sel(time='2017'), *normfit)

# fit gev
gevfit = genextreme.fit(qn_detrended_2017r.values)
qn_detrended_2017r_tau_parametric_gev = 1/genextreme.sf(qn.sel(time='2017'), *gevfit)

# fit gumbel_r
gumrfit = gumbel_r.fit(qn_detrended_2017r.values)
qn_detrended_2017r_tau_parametric_gum = 1/gumbel_r.sf(qn.sel(time='2017'), *gumrfit)

df.loc['QN 2017 parametric tau (Norm)','2017 removed linear detrended'] = qn_detrended_2017r_tau_parametric_norm.squeeze().round(2)
df.loc['QN 2017 parametric tau (GEV)','2017 removed linear detrended'] = qn_detrended_2017r_tau_parametric_gev.squeeze().round(2)
df.loc['QN 2017 parametric tau (Gumbel)','2017 removed linear detrended'] = qn_detrended_2017r_tau_parametric_gum.squeeze().round(2)

# get trend parameters
df.loc['QN 2017 predict', '2017 removed linear model'] = (dtrend_2017r['b']*2017+dtrend_2017r['a']).round(2)
df.loc['slope (ºC/year)', '2017 removed linear model'] = dtrend_2017r['b'].round(2)
df.loc['intercept (year)', '2017 removed linear model'] = dtrend_2017r['a'].round(2)
df.loc['R2', '2017 removed linear model'] = (dtrend_2017r['r']**2).round(2)

df.to_csv('../../../megafires_data/output/data.csv')


