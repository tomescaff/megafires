import xarray as xr
import numpy as np
from scipy import stats
import pandas as pd

def from_mon_to_year(series_mon_ave, init_date='1979-01-01', end_date='2006-12-31'):

    dr = pd.date_range(init_date, end_date, freq='1D')
    one_per_day = xr.DataArray(np.ones((dr.size,)), coords=[dr], dims=['time'])
    days_per_month = one_per_day.resample(time='1MS').sum(skipna=False)
    days_per_year = one_per_day.resample(time='1YS').sum(skipna=False)

    series_mon_acc = series_mon_ave*days_per_month
    series_yea_acc = series_mon_acc.resample(time='1YS').sum(skipna=False)
    series_yea_ave = series_yea_acc/days_per_year

    return series_yea_ave


ds = xr.open_dataset('../../../megafires_data/LENS_ALL/LENS_tmean_mean_mon.nc')
da = ds['TSA']-273.15
da_seltime = da.sel(time=slice('1979', '2006'))
da_year = from_mon_to_year(da_seltime)

ds_out = xr.Dataset({'TSA':da_year})
ds_out.to_netcdf('../../../megafires_data/LENS_ALL/LENS_tmean_1979_2006_yearly.nc')
