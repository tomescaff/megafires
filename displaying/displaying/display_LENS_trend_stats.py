import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import xarray as xr
import numpy as np
import cmaps
import sys

sys.path.append('../../processing')

import processing.utils as ut

# get LENS trend 1950-2020 jan tmax
ds = xr.open_dataset('../../../megafires_data/LENS_ALL/LENS_tmax_mean_mon_1950_2020_jan_trend.nc')
da = ds['slope'].mean('run')*10

# lat lon point of interests
qnll = (-70.6828, -33.4450) # quinta normal
cull = (-71.2167, -34.9664) # curico
chll = (-72.0400, -36.5872) # chillan


# create figure
fig = plt.figure(figsize=(8,7))

# define projection
ax = plt.axes(projection=ccrs.PlateCarree())

# set extent of map
ax.set_extent([-77, -66.5, -56.5, -17], crs=ccrs.PlateCarree())

# define and set  x and y ticks
xticks = [ -76, -72, -68]
yticks = [ -55, -50, -45, -40, -35, -30, -25, -20]
ax.set_xticks( xticks, crs=ccrs.PlateCarree())
ax.set_yticks( yticks, crs=ccrs.PlateCarree())

# format x and y labels
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

# add backgroung for land and ocean
resol = '50m'  
land = cfeature.NaturalEarthFeature('physical', 'land',  scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land'])
ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='none', facecolor=cfeature.COLORS['water'])

ax.add_feature(land, linewidth=0.0, alpha=0.5)
ax.add_feature(ocean, alpha = 0.5)

# add grid using previous ticks
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='grey', alpha=0.7, linestyle='--', draw_labels=False)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)

# plot the climatology and reshape color bar
pcm = ax.pcolormesh(da.lon.values, da.lat.values, da.values, cmap='jet', alpha=0.8, zorder=3, edgecolor='k', lw=0.1, vmin=0.0, vmax=0.5)
cbar = plt.colorbar(pcm, aspect = 40, pad=0.03)

lats = [-34.397907, -33.455498, -32.51309]
lons = [-71.25,  -70.0]
xx, yy = np.meshgrid(lons, lats)
ax.scatter( np.ravel(xx), np.ravel(yy), zorder = 4, edgecolor='r', facecolor='none')

ax.scatter( *qnll, edgecolor='grey', facecolor='none', zorder=4)
ax.scatter( *cull, edgecolor='grey', facecolor='none', zorder=4)
ax.scatter( *chll, edgecolor='grey', facecolor='none', zorder=4)

# draw the coastlines
land = cfeature.NaturalEarthFeature('physical', 'land',  scale=resol, edgecolor='k', facecolor='none')
ax.add_feature(land, linewidth=0.5, alpha=1, zorder=4)

#  reduce outline patch linewidths
cbar.outline.set_linewidth(0.4)
ax.outline_patch.set_linewidth(0.4)

# reduce fontsize
cbar.ax.tick_params(labelsize=8) 

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
    
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(8) 

# set title
cbar.ax.get_yaxis().labelpad = 12
cbar.ax.set_ylabel('LENS trend (ensmean) january Tmax 1950-2020 (ºC/dec)', fontdict={'fontsize':10})
plt.savefig('../../../megafires_data/png/LENS_trend_ensmean_jan_tmax_1950_2020_with_nearest_gridpoints.png', dpi=300)
plt.show()