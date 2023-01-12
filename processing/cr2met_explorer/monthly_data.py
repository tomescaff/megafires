import pandas as pd
import numpy as np
import xarray as xr

def daily_mean_with_nans(x, axis, **kwargs):
    frac=0.8
    cnt = lambda col: np.count_nonzero(~np.isnan(col))
    red = lambda col : np.nan if cnt(col)/col.size < frac else np.nanmean(col)
    ans = np.apply_along_axis(red, axis, x)
    return ans

def monthly_mean_with_nans(x):
    frac = 0.8
    cnt = lambda col: np.count_nonzero(~np.isnan(col))
    red = lambda col : np.nan if cnt(col)/col.size < frac else np.nanmean(col)
    ans = red(x)
    return ans

stnfile = '../../../megafires_data/stns/cr2_tasmaxDaily_2022.nc'

ds = xr.open_dataset(stnfile)
ds = ds.where( (ds.lat <= -26) & (ds.lat >= -40), drop=True)

tmax = ds['tmax']
lat = ds['lat']
lon = ds['lon']
alt = ds['alt']
name = ds['name']

tmax = tmax.where(tmax != -9999.0)
tmax_mon = tmax.resample(time='1MS').reduce(daily_mean_with_nans, dim='time')
tmax_jan = tmax_mon[::12]
tmax_jan = tmax_jan.sel(time=slice('1930', '2021'))

# for the anom map
clim_period = tmax_jan.sel(time=slice('1991', '2020'))
clim = xr.apply_ufunc(monthly_mean_with_nans, clim_period, input_core_dims=[["time"]], vectorize = True)
anom = tmax_jan.sel(time='2017') - clim
anom = anom.dropna('stn')

# fro the return period map
mat = np.zeros((tmax_jan.stn.size,))
for i in range(tmax_jan.stn.size):
    k = 0
    for j in range(tmax_jan.time.size):
        if np.isnan(tmax_jan[-(j+1),i]):
            break
        else:
            k=k+1
    mat[i] = k

mask = mat>31
tmax_jan_long_record = tmax_jan[:, mask]
lat_long_record = lat[mask]
lon_long_record = lon[mask]
alt_long_record = alt[mask]
name_long_record = name[mask]






