import xarray as xr
import numpy as np
from scipy import stats

# read LENS netcdf
ds = xr.open_dataset('../../../megafires_data/LENS2_ALL/LENS2_tmax_mean_mon.nc')
da = ds['TREFMXAV'] - 273.15
da = da.where(da.time.dt.month==1, drop=True)
da.to_netcdf('../../../megafires_data/LENS2_ALL/LENS2_tmax_mean_mon_jan.nc')

