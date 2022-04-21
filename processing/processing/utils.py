import pandas as pd
import xarray as xr
import numpy as np
from scipy import stats

# get Quinta Normal time series of tmax january monthly mean
def get_QN_series():
    
    # define path of csv file with data
    filepath = '../../../megafires_data/series/tmax_january_mon_mean_QN.csv'

    # read the csv file using pandas
    df = pd.read_csv(filepath, delimiter=",", decimal=",", parse_dates={'time': ['agno', 'mes', 'dia']})
    df = df.rename({'valor':'tmax'}, axis='columns')
    df = df.set_index('time')

    # to xarray data array
    da = df['tmax'].to_xarray()

    return da

# get CR2MET tmax over Chilean territory from 1979 to 2018 from file
def get_CR2MET():

    # define path of netcdf file with data
    filepath = '../../../megafires_data/CR2MET/CR2MET_tmax_v2.0_mon_1979_2018_005deg.nc'

    # open dataset
    ds = xr.open_dataset(filepath, decode_times=False)

    # reassign time
    ds['time'] = pd.date_range('1979-01', '2018-12', freq='1MS')

    return ds['tmax']

# get CR2MET january tmax 
def get_CR2MET_jan():

    da = get_CR2MET()

    return da[::12, :,:]

# get CR2MET corr jan with QN jan tmax
def get_corr_CR2MET_QN_jan_tmax():

    # get QN time series
    qn = get_QN_series()
    qn = qn.sel(time=slice('1979-01', '2018-12'))

    # get CR2MET jan tmax
    cr2 = get_CR2MET_jan()

    # to numpy array
    x = cr2.values
    y = qn.values

    t, n, m = x.shape

    # init empty matrix 
    matrix = np.zeros((n, m))
    
    # fill matrix with correlations
    for i in range(n):
        for j in range(m):
            if np.isnan(x[:,i,j]).all() or np.isnan(y[:]).all():
                matrix[i,j] = np.nan
            elif (x[:,i,j]==x[0,i,j]).all() or (y[:]==y[0]).all():
                matrix[i,j] = np.nan
            else:
                r,p = stats.pearsonr(x[:,i,j], y)
                matrix[i,j] = r
    
    # to xarray DataArray
    da = xr.DataArray(data = matrix, coords = [cr2.lat, cr2.lon], dims = ['lat', 'lon'])
    
    return da

# get CR2MET 2017 anomalies (clim. 1981-2010)
def get_CR2MET_2017_anom():

    # get CR2MET jan tmax
    cr2 = get_CR2MET_jan()

    # compute clim
    clim = cr2.sel(time=slice('1981-01-01', '2010-12-31')).mean('time').squeeze()
    
    # get 2017 value
    value = cr2.sel(time=slice('2017-01-01', '2017-12-31')).squeeze()

    # compute anom
    anom = value - clim 

    return anom

