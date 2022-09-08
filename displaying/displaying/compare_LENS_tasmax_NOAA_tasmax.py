import xarray as xr
import matplotlib.pyplot as plt

basedir = '/Volumes/Magister/tcarrasc/datasets/attribution/models/'

noaa = xr.open_dataset(basedir+'tasmax_CESM1-CAM5_LENS_mean_QN_latlon.nc')['tasmax'].squeeze()
lens = xr.open_dataset(basedir+'LENS_tasmax_mon_mean_QN_latlon.nc')['TREFMXAV']
lensmean = lens.mean('run').squeeze()


fig = plt.figure()
(noaa-lensmean).plot(ax=plt.gca(), color='r')
plt.show()

# no son iguales