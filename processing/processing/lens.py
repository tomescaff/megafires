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

def get_LENS_jan_tmax_QNW_control_run():
    # nearerst neighbor
    basedir = '../../../megafires_data/LENS/'
    filename = 'LENS_tasmax_mon_mean_QN_control_run.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']
    return da

def get_LENS_jan_tmax_QNWE_control_run():
    # FIX: this is the version for NN
    basedir = '../../../megafires_data/LENS/'
    filename = 'LENS_tasmax_mon_mean_QN_control_run.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']
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

