import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.lens as lens
import processing.math as pmath
import processing.stations as stns

lens2_gmst_full = gmst.get_gmst_annual_lens2_ensmean()
lens2_tmax_full = lens.get_LENS2_jan_tmax_QNWE()
qn_full = stns.get_QN_tmax_jan()

qn = qn_full.sel(time=slice('1930', '2021'))
lens2_gmst = lens2_gmst_full.sel(time=slice('1930', '2021'))
lens2_tmax = lens2_tmax_full.sel(time=slice('1930', '2021'))

lens2_tmax_anom = lens2_tmax - lens2_tmax.mean('time')
lens2_tmax_corrected = lens2_tmax_anom + qn.mean('time')

lens2_gmst_arr = np.tile(lens2_gmst.values, lens2_tmax_corrected.shape[0])
lens2_tmax_arr = np.ravel(lens2_tmax_corrected.values)

# # bootstrap
nboot = 1000#10000
bspreds_mu0 = np.zeros((nboot,))
bspreds_sigma0 = np.zeros((nboot,))
bspreds_alpha = np.zeros((nboot,))

for i in range(nboot):
    qn_i, sm_i = bootstrap(lens2_tmax_arr, lens2_gmst_arr)
    xopt_i = pmath.mle_norm_2d_fast(qn_i, sm_i, [29.57, 1.15, 0.88])
    bspreds_mu0[i] = xopt_i[0]
    bspreds_sigma0[i] = xopt_i[1]
    bspreds_alpha[i] = xopt_i[2]

iter = np.arange(nboot)
ds = xr.Dataset({
    'mu0':    xr.DataArray(bspreds_mu0,    coords=[iter], dims=['iter']),
    'sigma0': xr.DataArray(bspreds_sigma0, coords=[iter], dims=['iter']),
    'alpha':  xr.DataArray(bspreds_alpha,  coords=[iter], dims=['iter']), 
})
filepath = '../../../megafires_data/output/MLE_tasmax_jan_LENS2_GMST_'+str(nboot)+'_normal_validation.nc'
ds.to_netcdf(filepath)