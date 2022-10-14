import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.lens as lens
import processing.math as pmath

lens1_gmst_full = gmst.get_gmst_annual_lens1_ensmean()
lens1_tmax_full = lens.get_LENS_jan_tmax_QNW()

lens1_gmst = lens1_gmst_full.sel(time=slice('2010', '2100'))
lens1_tmax = lens1_tmax_full.sel(time=slice('2010', '2100'))

lens1_gmst_arr = np.tile(lens1_gmst.values, lens1_tmax.shape[0])
lens1_tmax_arr = np.ravel(lens1_tmax.values)

# # bootstrap
nboot = 1000
bspreds_mu0 = np.zeros((nboot,))
bspreds_sigma0 = np.zeros((nboot,))
bspreds_alpha = np.zeros((nboot,))

for i in range(nboot):
    qn_i, sm_i = bootstrap(lens1_tmax_arr, lens1_gmst_arr)
    xopt_i = pmath.mle_norm_2d_fast(qn_i, sm_i, [29.55, 1.03, 1.11])
    bspreds_mu0[i] = xopt_i[0]
    bspreds_sigma0[i] = xopt_i[1]
    bspreds_alpha[i] = xopt_i[2]

iter = np.arange(nboot)
ds = xr.Dataset({
    'mu0':    xr.DataArray(bspreds_mu0,    coords=[iter], dims=['iter']),
    'sigma0': xr.DataArray(bspreds_sigma0, coords=[iter], dims=['iter']),
    'alpha':  xr.DataArray(bspreds_alpha,  coords=[iter], dims=['iter']), 
})
filepath = '../../../megafires_data/output/MLE_tasmax_jan_LENS1_GMST_'+str(nboot)+'_normal_future_QN_NN.nc'
ds.to_netcdf(filepath)