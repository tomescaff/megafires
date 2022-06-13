import xarray as xr
import numpy as np
import pandas as pd

# get Quinta Normal time series of tmax january monthly mean
def get_QN_tmax_jan_from_CR2():
    
    # define path of csv file with data
    filepath = '../../../megafires_data/series/tmax_january_mon_mean_QN.csv'

    # read the csv file using pandas
    df = pd.read_csv(filepath, delimiter=",", decimal=",", parse_dates={'time': ['agno', 'mes', 'dia']})
    df = df.rename({'valor':'tmax'}, axis='columns')
    df = df.set_index('time')

    # to xarray data array
    da = df['tmax'].to_xarray()

    return da

# get Quinta Normal time series of tmax january monthly mean from 1911 to 2013
def get_QN_tmax_jan_from_paceituno():

    # define path of csv file with data
    filepath = '../../../megafires_data/QN/T_maxima_diaria_QN_1911_2013.csv'

    # read the csv file using pandas
    df = pd.read_csv(filepath, decimal=",")

    # select january tmax
    df_jan = df.iloc[0:31,:]

    # create xarray with mean values rounded
    years = np.arange(1911,2014)
    data = np.zeros((years.size,))

    for i, year in enumerate(years):
        data[i] = round(df_jan.loc[:,f'{year}'].mean(), 1)

    time = pd.date_range('1911-01-01', '2013-01-01',freq='1YS')

    da = xr.DataArray(data, coords=[time], dims=['time'])
    return da

# get Quinta Normal time series of tmax january monthly mean from 1911 to 2021
def get_QN_tmax_jan():

    cr2 = get_QN_tmax_jan_from_CR2()
    paceituno = get_QN_tmax_jan_from_paceituno()
    paceituno_slice = paceituno.sel(time=slice('1911', '1949'))

    full = xr.concat([paceituno_slice, cr2], dim='time')
    return full

# get Curico time series of tmax january monthly mean
def get_CU_tmax_jan():
    
    # define path of csv file with data
    filepath = '../../../megafires_data/series/tmax_january_mon_mean_CUDMC.csv'

    # read the csv file using pandas
    df = pd.read_csv(filepath, delimiter=",", decimal=".", parse_dates={'time': ['agno', ' mes', ' dia']})
    df = df.rename({' valor':'tmax'}, axis='columns')
    df = df.set_index('time')

    # to xarray data array
    da = df['tmax'].to_xarray()

    return da

# get Chillan time series of tmax january monthly mean
def get_CH_tmax_jan():
    
    # define path of csv file with data
    filepath = '../../../megafires_data/series/tmax_january_mon_mean_CHDMC.csv'

    # read the csv file using pandas
    df = pd.read_csv(filepath, delimiter=",", decimal=".", parse_dates={'time': ['agno', ' mes', ' dia']})
    df = df.rename({' valor':'tmax'}, axis='columns')
    df = df.set_index('time')

    # to xarray data array
    da = df['tmax'].to_xarray()

    return da

def get_regindex_tmax_jan(norm=False, anom=True):

    qn = get_QN_tmax_jan().sel(time=slice('1966', '2021'))
    cu = get_CU_tmax_jan().sel(time=slice('1966', '2021'))
    ch = get_CH_tmax_jan().sel(time=slice('1966', '2021'))

    qn_1981_2010_mean = qn.sel(time=slice('1981', '2010')).mean('time')
    cu_1981_2010_mean = cu.sel(time=slice('1981', '2010')).mean('time')
    ch_1981_2010_mean = ch.sel(time=slice('1981', '2010')).mean('time')

    qn_1981_2010_std = qn.sel(time=slice('1981', '2010')).std('time')
    cu_1981_2010_std = cu.sel(time=slice('1981', '2010')).std('time')
    ch_1981_2010_std = ch.sel(time=slice('1981', '2010')).std('time')

    qn_anom = qn - qn_1981_2010_mean
    cu_anom = cu - cu_1981_2010_mean
    ch_anom = ch - ch_1981_2010_mean

    regindex_anom = (qn_anom+cu_anom+ch_anom)/3.0
    regindex_mean = (qn_1981_2010_mean + cu_1981_2010_mean + ch_1981_2010_mean)/3.0

    regindex_norm = (qn_anom/qn_1981_2010_std + cu_anom/cu_1981_2010_std + ch_anom/ch_1981_2010_std)/3.0

    if not norm:
        if anom:
            ans = regindex_anom
        else:
            ans = regindex_anom + regindex_mean
    else:
        ans = regindex_norm
    
    return ans


        

    




