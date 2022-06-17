import xarray as xr
import numpy as np
import pandas as pd

def get_burned_area_january():
    basedir = '../../../megafires_data/fires/'
    filename = 'dano_incendios_forestales_segun_mes.csv'
    filepath = basedir + filename
    df = pd.read_csv(filepath)
    ba_jan = df.iloc[7,1:38].astype(float)
    time = pd.date_range('1985-01', '2021-01', freq='1YS')
    da = xr.DataArray(ba_jan.values, coords=[time], dims=['time'])
    return da

def get_burned_area_february():
    basedir = '../../../megafires_data/fires/'
    filename = 'dano_incendios_forestales_segun_mes.csv'
    filepath = basedir + filename
    df = pd.read_csv(filepath)
    ba_feb = df.iloc[8,1:38].astype(float)
    time = pd.date_range('1985-02', '2021-02', freq='12MS')
    da = xr.DataArray(ba_feb.values, coords=[time], dims=['time'])
    return da

def get_burned_area_december():
    basedir = '../../../megafires_data/fires/'
    filename = 'dano_incendios_forestales_segun_mes.csv'
    filepath = basedir + filename
    df = pd.read_csv(filepath)
    ba_dec = df.iloc[6,1:38].astype(float)
    time = pd.date_range('1984-12', '2020-12', freq='12MS')
    da = xr.DataArray(ba_dec.values, coords=[time], dims=['time'])
    return da

def get_burned_area_by_season():
    basedir = '../../../megafires_data/fires/'
    filename = 'resumen_nacional_ocurrencia_dano_1964_2021.csv'
    filepath = basedir + filename
    df = pd.read_csv(filepath)
    burned_area = df.iloc[2:,4]
    time = pd.date_range('1964-01', '2021-01', freq='12MS')
    da = xr.DataArray(burned_area.values, coords=[time], dims=['time'])
    return da

def get_fires_frequency_by_season():
    basedir = '../../../megafires_data/fires/'
    filename = 'resumen_nacional_ocurrencia_dano_1964_2021.csv'
    filepath = basedir + filename
    df = pd.read_csv(filepath)
    fires_frequency = df.iloc[2:,3]
    time = pd.date_range('1964-01', '2021-01', freq='12MS')
    da = xr.DataArray(fires_frequency.values, coords=[time], dims=['time'])
    return da


 