import xarray as xr
import numpy as np

ds = xr.open_dataset('../../../megafires_data/LENS_ALL/LENS_tasmax_mon_mean_control_run.nc')
da = ds['TREFMXAV']
qne = da.sel(lat= -33.455497, lon = (-70.0)%360, method='nearest').drop(['lat','lon'])
qnw = da.sel(lat= -33.455497, lon = (-71.25)%360, method='nearest').drop(['lat','lon'])
qn = (qnw + qne)/2.0
qn = qn.where(qn.time.dt.month == 1, drop=True)
qn = qn - 273.15
dsout = xr.Dataset({'TREFMXAV': qn})
dsout.to_netcdf('../../../megafires_data/LENS/LENS_tasmax_mon_mean_QN_control_run.nc')