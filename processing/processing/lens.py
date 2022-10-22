import xarray as xr
import numpy as np

def get_LENS_jan_tmax_QNW():
    # nearest neighbor
    basedir = '../../../megafires_data/LENS/'
    filename = 'LENS_tmax_mon_QN.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']-273.15
    return da[:, ::12]

def get_LENS_jan_tmax_QNE():
    basedir = '../../../megafires_data/LENS/'
    filename = 'LENS_tmax_mon_QNE.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']-273.15
    return da[:, ::12]

def get_LENS_jan_tmax_QNWE():
    qnw = get_LENS_jan_tmax_QNW().drop(['lat','lon'])
    qne = get_LENS_jan_tmax_QNE().drop(['lat','lon'])
    return (qnw + qne)/2.0

def get_LENS_jan_tmax_CU_NN():
    # nearest neighbor
    basedir = '../../../megafires_data/LENS_ALL/'
    filename = 'LENS_tmax_mean_mon_chile.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']#-273.15
    lat, lon = -34.9664, -71.2167 
    da = da.sel(lat=lat, lon=lon%360, method = 'nearest')
    da = da.where(da.time.dt.month == 1, drop=True)
    return da

def get_LENS_jan_tmax_CH_NN():
    # nearest neighbor
    basedir = '../../../megafires_data/LENS_ALL/'
    filename = 'LENS_tmax_mean_mon_chile.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']#-273.15
    lat, lon = -36.5872, -72.0400
    da = da.sel(lat=lat, lon=lon%360, method = 'nearest')
    da = da.where(da.time.dt.month == 1, drop=True)
    return da

def get_LENS_jan_tmax_QNW_control_run():
    # nearerst neighbor
    basedir = '../../../megafires_data/LENS_ALL/'
    filename = 'LENS_tasmax_mon_mean_control_run_chile.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV'] - 273.15
    lon, lat = -70.6828, -33.4450
    da = da.sel(lat=lat, lon=lon%360, method = 'nearest')
    da = da.where(da.time.dt.month == 1, drop=True)
    return da

def get_LENS_jan_tmax_QNWE_control_run():
    # FIX: this is the version for NN
    basedir = '../../../megafires_data/LENS/'
    filename = 'LENS_tasmax_mon_mean_QN_control_run.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']
    return da

def get_LENS_jan_tmax_CU_control_run():
    # nearerst neighbor
    basedir = '../../../megafires_data/LENS_ALL/'
    filename = 'LENS_tasmax_mon_mean_control_run_chile.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV'] - 273.15
    lat, lon = -34.9664, -71.2167
    da = da.sel(lat=lat, lon=lon%360, method = 'nearest')
    da = da.where(da.time.dt.month == 1, drop=True)
    return da

def get_LENS_jan_tmax_CH_control_run():
    # nearerst neighbor
    basedir = '../../../megafires_data/LENS_ALL/'
    filename = 'LENS_tasmax_mon_mean_control_run_chile.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV'] - 273.15
    lat, lon = -36.5872, -72.0400
    da = da.sel(lat=lat, lon=lon%360, method = 'nearest')
    da = da.where(da.time.dt.month == 1, drop=True)
    return da

def get_LENS2_jan_tmax_QNW():
    # nearest neighbor
    basedir = '../../../megafires_data/LENS2/'
    filename = 'LENS2_tmax_mean_mon_QNW.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']-273.15
    return da[:, ::12]

def get_LENS2_jan_tmax_QNE():
    basedir = '../../../megafires_data/LENS2/'
    filename = 'LENS2_tmax_mean_mon_QNE.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']-273.15
    return da[:, ::12]

def get_LENS2_jan_tmax_QNWE():
    qnw = get_LENS2_jan_tmax_QNW().drop(['lat','lon'])
    qne = get_LENS2_jan_tmax_QNE().drop(['lat','lon'])
    return (qnw + qne)/2.0

def get_LENS2_jan_tmax_CU_NN():
    # nearest neighbor
    basedir = '../../../megafires_data/LENS2_ALL/'
    filename = 'LENS2_tmax_mean_mon.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']-273.15
    lat, lon = -34.9664, -71.2167 
    da = da.sel(lat=lat, lon=lon%360, method = 'nearest')
    da = da.where(da.time.dt.month == 1, drop=True)
    return da

def get_LENS2_jan_tmax_CH_NN():
    # nearest neighbor
    basedir = '../../../megafires_data/LENS2_ALL/'
    filename = 'LENS2_tmax_mean_mon.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']-273.15
    lat, lon = -36.5872, -72.0400
    da = da.sel(lat=lat, lon=lon%360, method = 'nearest')
    da = da.where(da.time.dt.month == 1, drop=True)
    return da

