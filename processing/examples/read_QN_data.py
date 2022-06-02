import pandas as pd
import numpy as np
import xarray as xr

df = pd.read_csv('../../../megafires_data/QN/T_maxima_diaria_QN_1911_2013.csv', decimal=",")
df_jan = df.iloc[0:31,:]

years = np.arange(1911,2014)
data = np.zeros((years.size,))

for i, year in enumerate(years):
    data[i] = round(df_jan.loc[:,f'{year}'].mean(), 1)

time = pd.date_range('1911-01-01', '2013-01-01',freq='1YS')

da = xr.DataArray(data, coords=[time], dims=['time'])