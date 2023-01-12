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

stnfile = '../../../megafires_data/CR2_explorer/cr2_tasmaxDaily_2022.nc'

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
anom_data = anom.squeeze()
anom_lats = lat.sel(stn=anom_data.stn)
anom_lons = lon.sel(stn=anom_data.stn)
anom_alts = alt.sel(stn=anom_data.stn)
anom_name = name.sel(stn=anom_data.stn)

ds_anom = xr.Dataset({'data':anom_data, 'lat':anom_lats, 'lon':anom_lons, 'alt':anom_alts, 'name':anom_name})
ds_anom.to_netcdf('../../../megafires_data/CR2_explorer/cr2_tasmax_anom2017_refperiod_1991_2020_26deg_40degS.nc', encoding={"name": {"dtype": "str"}})

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
long_record_data = tmax_jan[:, mask]
long_record_lats = lat[mask]
long_record_lons = lon[mask]
long_record_alts = alt[mask]
long_record_name = name[mask]

ds_long_record = xr.Dataset({'data':long_record_data, 'lat':long_record_lats, 'lon':long_record_lons, 'alt':long_record_alts, 'name':long_record_name})
ds_long_record.to_netcdf('../../../megafires_data/CR2_explorer/cr2_tasmax_stns_data_before_1990_26deg_40degS.nc', encoding={"name": {"dtype": "str"}})






