import xarray as xr
import numpy as np
from scipy import stats
import pandas as pd


ds = xr.open_dataset('../../../megafires_data/LENS_ALL/LENS_tmax_1979_2006_yearly.nc')
da = ds['TREFMXAV']

# get shape
nrun, ntime, nlat, nlon = da.values.shape

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
            series = da[ir,:,ilat,ilon]

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

da_slope = xr.DataArray(matrix_slope, coords = [da.run, da.lat, da.lon], dims=['run', 'lat', 'lon'])
da_intercept = xr.DataArray(matrix_intercept, coords = [da.run, da.lat, da.lon], dims=['run', 'lat', 'lon'])
da_rvalue = xr.DataArray(matrix_rvalue, coords = [da.run, da.lat, da.lon], dims=['run', 'lat', 'lon'])
da_pvalue = xr.DataArray(matrix_pvalue, coords = [da.run, da.lat, da.lon], dims=['run', 'lat', 'lon'])

ds_out = xr.Dataset({'slope':da_slope, 'intercept':da_intercept, 'rvalue':da_rvalue, 'pvalue':da_pvalue})
ds_out.to_netcdf('../../../megafires_data/LENS_ALL/LENS_tmax_1979_2006_yearly_trend.nc')
