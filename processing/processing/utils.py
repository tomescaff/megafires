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

# get Quinta Normal time series of tmax january monthly mean from 1911 to 2013
def get_QN_past_series():

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
def get_QN_full_series():

    cr2 = get_QN_series()
    paceituno = get_QN_past_series()
    paceituno_slice = paceituno.sel(time=slice('1911', '1949'))

    full = xr.concat([paceituno_slice, cr2], dim='time')
    return full

# get Curico time series of tmax january monthly mean
def get_CU_series():
    
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
def get_CH_series():
    
    # define path of csv file with data
    filepath = '../../../megafires_data/series/tmax_january_mon_mean_CHDMC.csv'

    # read the csv file using pandas
    df = pd.read_csv(filepath, delimiter=",", decimal=".", parse_dates={'time': ['agno', ' mes', ' dia']})
    df = df.rename({' valor':'tmax'}, axis='columns')
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

# get CR2MET corr jan with time series jan tmax
def get_corr_CR2MET_jan_tmax(series):

    # slice time series
    series = series.sel(time=slice('1979-01', '2018-12'))

    # get CR2MET jan tmax
    cr2 = get_CR2MET_jan()

    # to numpy array
    x = cr2.values
    y = series.values

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

# get CR2MET corr jan with QN jan tmax
def get_corr_CR2MET_QN_jan_tmax():

    qn = get_QN_series()
    return get_corr_CR2MET_jan_tmax(qn)
    
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

def get_CR2MET_2017_anom_std():
    
    # get CR2MET jan tmax
    cr2 = get_CR2MET_jan()

    # compute clim
    clim = cr2.sel(time=slice('1981-01-01', '2010-12-31')).mean('time').squeeze()

    # compute std
    std = cr2.sel(time=slice('1981-01-01', '2010-12-31')).std('time').squeeze()
    
    # get 2017 value
    value = cr2.sel(time=slice('2017-01-01', '2017-12-31')).squeeze()

    # compute anom
    anom = (value - clim)/std

    return anom

def get_CR2MET_2017_ranking():
    
    # get CR2MET jan tmax
    cr2 = get_CR2MET_jan()

    ntime, nlat, nlon = cr2.values.shape
    matrix = np.zeros((nlat,nlon))

    for ilat in range(nlat):
        for ilon in range(nlon):
            series = cr2[:, ilat, ilon].values
            if np.isnan(series[0]):
                matrix[ilat, ilon] = np.nan
            else:
                ranks = stats.rankdata(series)
                matrix[ilat, ilon] = ranks[-2]
    
    # to xarray DataArray
    da = xr.DataArray(data = matrix, coords = [cr2.lat, cr2.lon], dims = ['lat', 'lon'])
    
    return da

# get return period function from QN station data
def get_QN_tau():

    # get QN series
    da = get_QN_series()

    n = da.size

    # sort values
    z = da.values
    z = np.sort(z)

    # get unique values
    u = np.unique(z)

    m = u.size

    # create matrix for tail probability and tau
    tail = np.zeros((m,))
    tau = np.zeros((m,))

    # compute tail and tau
    for i in range(m):
        x = u[i]
        tail[i] = np.sum(z>=x)/n
        tau[i] = 1/tail[i]

    return u, tau

# get return period function from QN station data without 2017 max
def get_QN_tau_remove_max():

    # get QN series
    da = get_QN_series()

    n = da.size

    # sort values
    z = da.values
    z = np.sort(z)
    z = z[:-1]

    # get unique values
    u = np.unique(z)

    m = u.size

    # create matrix for tail probability and tau
    tail = np.zeros((m,))
    tau = np.zeros((m,))

    # compute tail and tau
    for i in range(m):
        x = u[i]
        tail[i] = np.sum(z>=x)/n
        tau[i] = 1/tail[i]

    return u, tau

# get return period function from parametric distribution
def get_param_tau(param_fun, args):

    # get random values
    z = param_fun(*args, size=1000) 

    n = z.size

    # sort values
    z = np.sort(z)

    # get unique values
    u = np.unique(z)

    m = u.size

    # create matrix for tail probability and tau
    tail = np.zeros((m,))
    tau = np.zeros((m,))

    # compute tail and tau
    for i in range(m):
        x = u[i]
        tail[i] = np.sum(z>=x)/n
        tau[i] = 1/tail[i]

    return u, tau

# get linear trend from QN station data
def get_linear_trend():

    # get QN series
    da = get_QN_series()
    y = da.values
    x = da.time.dt.year.values

    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

    y_pred = slope*x + intercept

    ans = dict()
    ans['y'] = y
    ans['x'] = x
    ans['y_pred'] = y_pred
    ans['b'] = slope
    ans['a'] = intercept
    ans['r'] = r_value
    ans['p'] = p_value
    ans['std_err'] = std_err

    return ans

# get linear trend from QN station data with 2017 removed
def get_linear_trend_2017r():

    # get QN series
    da = get_QN_series()
    da = da.where(da.time.dt.year != 2017, drop = True)
    y = da.values
    x = da.time.dt.year.values

    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

    y_pred = slope*x + intercept

    ans = dict()
    ans['y'] = y
    ans['x'] = x
    ans['y_pred'] = y_pred
    ans['b'] = slope
    ans['a'] = intercept
    ans['r'] = r_value
    ans['p'] = p_value
    ans['std_err'] = std_err

    return ans

# get return period function from QN station data detrended
def get_QN_tau_detrended():

    # get QN series
    da = get_QN_series()

    # get linear trend parameters from full series
    dtrend = get_linear_trend()

    # get predicted values from full series
    qn_pred = xr.DataArray(dtrend['y_pred'], coords=[da.time], dims=['time'])

    # get detrended series
    qn_detrended = da - qn_pred + da.mean('time')

    n = qn_detrended.size

    # sort values
    z = qn_detrended.values
    z = np.sort(z)

    # get unique values
    u = np.unique(z)

    m = u.size

    # create matrix for tail probability and tau
    tail = np.zeros((m,))
    tau = np.zeros((m,))

    # compute tail and tau
    for i in range(m):
        x = u[i]
        tail[i] = np.sum(z>=x)/n
        tau[i] = 1/tail[i]

    return u, tau

# get Quinta Normal time series of tmax january monthly mean detrended
def get_QN_series_detrended():

    # get QN series
    da = get_QN_series()

    # get linear trend parameters from full series
    dtrend = get_linear_trend()

    # get predicted values from full series
    qn_pred = xr.DataArray(dtrend['y_pred'], coords=[da.time], dims=['time'])

    # get detrended series
    qn_detrended = da - qn_pred + da.mean('time')

    return qn_detrended

# get Quinta Normal time series of tmax january monthly mean detrended 2017 removed
def get_QN_series_detrended_2017r():

    # get QN series
    da = get_QN_series()

    # remove 2017
    da_2017r = da.where(da.time.dt.year != 2017, drop = True)

    # get linear trend parameters from 2017 removed series
    dtrend_2017r = get_linear_trend_2017r()

    # get predicted values from full series
    da_pred_2017r = xr.DataArray(dtrend_2017r['y_pred'], coords=[da_2017r.time], dims=['time'])

    # get detrended series
    da_detrended_2017r = da_2017r - da_pred_2017r + da_2017r.mean('time')

    return da_detrended_2017r

# get return period function from QN station data detrended 2017 removed
def get_QN_tau_detrended_2017r():

    # get detrended series
    qn_detrended_2017r = get_QN_series_detrended_2017r()

    n = qn_detrended_2017r.size

    # sort values
    z = qn_detrended_2017r.values
    z = np.sort(z)

    # get unique values
    u = np.unique(z)

    m = u.size

    # create matrix for tail probability and tau
    tail = np.zeros((m,))
    tau = np.zeros((m,))

    # compute tail and tau
    for i in range(m):
        x = u[i]
        tail[i] = np.sum(z>=x)/n
        tau[i] = 1/tail[i]

    return u, tau

# get LENS jan Tmax time series
def get_LENS_jan():
    # get LENS time series
    basedir = '../../../megafires_data/LENS/'
    filename = 'LENS_tmax_mon_QN.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']-273.15
    da_jan = da[:, ::12]
    return da_jan
    
# get LENS jan Tmax time series 40 runs as one numpy array between 1950-2021
def get_LENS_jan_1950_2021_ravel():
    da_jan = get_LENS_jan()
    da_jan_1950_2021 = da_jan.sel(time=slice('1950', '2021'))
    np_jan_all_runs_1950_2021 = np.ravel(da_jan_1950_2021.values)
    return np_jan_all_runs_1950_2021

# get control run jan tmax as numpy array
def get_control_run_jan():
    basedir = '../../../megafires_data/LENS/'
    filename = 'LENS_tmax_mon_QN_control_run.nc'
    filepath = basedir + filename
    ds = xr.open_dataset(filepath)
    da = ds['TREFMXAV']-273.15
    da_jan = da[::12]
    np_jan_control_run = np.ravel(da_jan.values)
    return np_jan_control_run

# get return periods from arg as numpy array
def get_LENS_jan_1950_2021_tau(z):
    
    n = z.size

    # sort values
    z = np.sort(z)

    # get unique values
    u = np.unique(z)

    m = u.size

    # create matrix for tail probability and tau
    tail = np.zeros((m,))
    tau = np.zeros((m,))

    # compute tail and tau
    for i in range(m):
        x = u[i]
        tail[i] = np.sum(z>=x)/n
        tau[i] = 1/tail[i]

    return u, tau

def get_tau(z):
    
    n = z.size

    # sort values
    z = np.sort(z)

    # get unique values
    u = np.unique(z)

    m = u.size

    # create matrix for tail probability and tau
    tail = np.zeros((m,))
    tau = np.zeros((m,))

    # compute tail and tau
    for i in range(m):
        x = u[i]
        tail[i] = np.sum(z>=x)/n
        tau[i] = 1/tail[i]

    return u, tau

# compute annual trend from [run, time, lat, lon] xarray DataArray
def compute_trend(xr_array):

    # get shape
    nrun, ntime, nlat, nlon = xr_array.values.shape

    # initalize matrix with results
    matrix_slope = np.zeros((nrun,nlat,nlon))
    matrix_intercept = np.zeros((nrun,nlat,nlon))
    matrix_rvalue = np.zeros((nrun,nlat,nlon))
    matrix_pvalue = np.zeros((nrun,nlat,nlon))

    # for each run-lat-lon
    for ir in range(nrun):
        for ilat in range(nlat):
            for ilon in range(nlon):

                # get the series
                series = xr_array[ir,:,ilat,ilon]

                # get the x and y values
                y = series.values
                x = series.time.dt.year.values

                # compute the linear regression of y onto x
                slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

                # save results
                matrix_slope[ir, ilat, ilon] = slope
                matrix_intercept[ir, ilat, ilon] = intercept
                matrix_rvalue[ir, ilat, ilon] = r_value
                matrix_pvalue[ir, ilat, ilon] = p_value

    da_slope = xr.DataArray(matrix_slope, coords = [xr_array.run, xr_array.lat, xr_array.lon], dims=['run', 'lat', 'lon'])
    da_intercept = xr.DataArray(matrix_intercept, coords = [xr_array.run, xr_array.lat, xr_array.lon], dims=['run', 'lat', 'lon'])
    da_rvalue = xr.DataArray(matrix_rvalue, coords = [xr_array.run, xr_array.lat, xr_array.lon], dims=['run', 'lat', 'lon'])
    da_pvalue = xr.DataArray(matrix_pvalue, coords = [xr_array.run, xr_array.lat, xr_array.lon], dims=['run', 'lat', 'lon'])
    ds_out = xr.Dataset({'slope':da_slope, 'intercept':da_intercept, 'rvalue':da_rvalue, 'pvalue':da_pvalue})
    return ds_out