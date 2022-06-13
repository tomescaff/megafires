import xarray as xr
import numpy as np
from scipy import stats

# read LENS netcdf
ds = xr.open_dataset('../../../megafires_data/LENS2_ALL/LENS2_tmax_mean_mon.nc')
ds_eastern = ds.sel(lat= -33.455498, lon = (-70.0)%360)
ds_eastern.to_netcdf('../../../megafires_data/LENS2/LENS2_tmax_mean_mon_QNE.nc')

ds_western = ds.sel(lat= -33.455498, lon = (-71.25)%360)
ds_western.to_netcdf('../../../megafires_data/LENS2/LENS2_tmax_mean_mon_QNW.nc')

