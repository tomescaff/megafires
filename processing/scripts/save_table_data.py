import sys
import pandas as pd

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

# compute taus with 2017 removed

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

df.to_csv('../../../megafires_data/output/data.csv')

