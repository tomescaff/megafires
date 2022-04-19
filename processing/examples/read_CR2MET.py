import xarray as xr
import pandas as pd

# define path of netcdf file with data
filepath = '../../../megafires_data/CR2MET/CR2MET_tmax_v2.0_mon_1979_2018_005deg.nc'

# open dataset
ds = xr.open_dataset(filepath, decode_times=False)

# reassign time
ds['time'] = pd.date_range('1979-01', '2018-12', freq='1MS')
