import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from sklearn.utils import resample as bootstrap

sys.path.append('../../processing')

import processing.gmst as gmst
import processing.math as pmath


smf = gmst.get_gmst_annual_5year_smooth()
stnsf = xr.open_dataset('../../../megafires_data/CR2_explorer/cr2_tasmax_stns_data_before_1990_26deg_40degS.nc')

N = stnsf.stn.size

for j in range(N):
    tmaxf = stnsf.data[:, j]
    name = str(stnsf.name[j].values)

    prev_year = int(xr.where(np.isnan(tmaxf), 1,np.nan).dropna('time').time.dt.year[-1].values)
    tx = tmaxf.sel(time=slice(str(prev_year+1), '2021'))
    sm = smf.sel(time=slice(str(prev_year+1), '2021'))
    
    sm = sm.where(sm.time.dt.year != 2017, drop=True)
    tx = tx.where(tx.time.dt.year != 2017, drop=True)

    # # bootstrap
    nboot = 10
    bspreds_mu0 = np.zeros((nboot,))
    bspreds_sigma0 = np.zeros((nboot,))
    bspreds_alpha = np.zeros((nboot,))

    for i in range(nboot):
        tx_i, sm_i = bootstrap(tx.values, sm.values)
        xopt_i = pmath.mle_norm_2d(tx_i, sm_i, [29.44, 0.81, 1.41])
        bspreds_mu0[i] = xopt_i[0]
        bspreds_sigma0[i] = xopt_i[1]
        bspreds_alpha[i] = xopt_i[2]

    iter = np.arange(nboot)
    ds = xr.Dataset({
        'mu0':    xr.DataArray(bspreds_mu0,    coords=[iter], dims=['iter']),
        'sigma0': xr.DataArray(bspreds_sigma0, coords=[iter], dims=['iter']),
        'alpha':  xr.DataArray(bspreds_alpha,  coords=[iter], dims=['iter']), 
    })
    filepath = f'../../../megafires_data/output/MLE_tasmax_jan_GMST_{nboot}_normal_validation_stn_{j:02}.nc'
    ds.to_netcdf(filepath)