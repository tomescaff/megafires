import pandas as pd
import xarray as xr

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
