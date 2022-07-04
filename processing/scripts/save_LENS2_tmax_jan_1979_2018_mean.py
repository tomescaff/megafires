import xarray as xr
import numpy as np
from scipy import stats

# read LENS netcdf
ds = xr.open_dataset('../../../megafires_data/LENS2_ALL/LENS2_tmax_mean_mon_jan.nc')
da = ds['TREFMXAV'].sel(time=slice('1979', '2018'))
da_mean = da.mean(['time','run']).squeeze()
dsout = xr.Dataset({'TREFMXAV':da_mean})
dsout.to_netcdf('../../../megafires_data/LENS2_ALL/LENS2_tmax_mean_mon_jan_1979_2018_mean.nc')

