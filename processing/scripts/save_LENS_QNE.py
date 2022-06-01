import xarray as xr
import numpy as np
from scipy import stats

# read LENS netcdf
ds = xr.open_dataset('../../../megafires_data/LENS_ALL/LENS_tmax_mean_mon.nc')
ds = ds.sel(lat= -33.455498, lon = (-70.0)%360)
ds.to_netcdf('../../../megafires_data/LENS/LENS_tmax_mean_mon_QNE.nc')

