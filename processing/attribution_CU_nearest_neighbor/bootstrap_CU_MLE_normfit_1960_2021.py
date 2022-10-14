import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.stations as stns
import processing.math as pmath


smf = gmst.get_gmst_annual_5year_smooth()
qnf = stns.get_CU_tmax_jan()

sm = smf.sel(time=slice('1960','2021'))
qn = qnf.sel(time=slice('1960','2021'))

sm = sm.where(sm.time.dt.year != 2017, drop=True)
qn = qn.where(qn.time.dt.year != 2017, drop=True)

# # bootstrap
nboot = 1000
bspreds_mu0 = np.zeros((nboot,))
bspreds_sigma0 = np.zeros((nboot,))
bspreds_alpha = np.zeros((nboot,))

for i in range(nboot):
    qn_i, sm_i = bootstrap(qn.values, sm.values)
    xopt_i = pmath.mle_norm_2d(qn_i, sm_i, [29.44, 0.81, 1.41])
    bspreds_mu0[i] = xopt_i[0]
    bspreds_sigma0[i] = xopt_i[1]
    bspreds_alpha[i] = xopt_i[2]

iter = np.arange(nboot)
ds = xr.Dataset({
    'mu0':    xr.DataArray(bspreds_mu0,    coords=[iter], dims=['iter']),
    'sigma0': xr.DataArray(bspreds_sigma0, coords=[iter], dims=['iter']),
    'alpha':  xr.DataArray(bspreds_alpha,  coords=[iter], dims=['iter']), 
})
filepath = '../../../megafires_data/output/MLE_tasmax_jan_CU_GMST_'+str(nboot)+'_normal_validation.nc'
ds.to_netcdf(filepath)