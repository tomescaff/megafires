import xarray as xr

# read LENS netcdf
ds = xr.open_dataset('../../../megafires_data/LENS_ALL/LENS_tmax_mean_mon.nc')
da = ds['TREFMXAV'] - 273.15
da = da.sel(lat=slice(-56.5,-17), lon=slice(283,293.5))
dsout = xr.Dataset({'TREFMXAV':da})
dsout.to_netcdf('../../../megafires_data/LENS_ALL/LENS_tmax_mean_mon_chile.nc')