import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io.shapereader import Reader

import sys
import cmaps
from scipy.stats import linregress
import numpy as np
import xarray as xr

sys.path.append('../../processing')

import processing.utils as ut

# get CR2MET january tmax series
cr2 = ut.get_CR2MET_jan()

# compute trend
da = cr2.sel(time=slice('1979-01-01', '2018-12-31'))

t, n, m = da.values.shape
matrix = np.zeros((n,m))
for i in range(n):
    for j in range(m):
        series = da[:,i,j]
        if np.isnan(series.values[0]):
            matrix[i,j] = np.nan
        else:
            slope, intercept, rval, pval, stderr = linregress(series.time.dt.year, series.values)
            matrix[i,j] = slope

da_trend = xr.DataArray(matrix, coords=[da.lat, da.lon], dims=['lat','lon'])

fname = '../../../megafires_data/shp/Regiones/Regional.shp'

# create figure
fig = plt.figure(figsize=(8,7))

# define projection
ax = plt.axes(projection=ccrs.PlateCarree())

# set extent of map
ax.set_extent([-76, -66.5, -56.5, -17], crs=ccrs.PlateCarree())

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

# plot the climatology and reshape color bar cmaps.amwg_blueyellowred
pcm = ax.pcolormesh(da_trend.lon.values, da_trend.lat.values, da_trend.values*100, cmap=cmaps.amwg_blueyellowred, zorder=4, vmin=-6, vmax=6)
cbar = plt.colorbar(pcm, aspect = 40, pad=0.03)

# draw the coastlines
land = cfeature.NaturalEarthFeature('physical', 'land',  scale=resol, edgecolor='k', facecolor='none')
ax.add_feature(land, linewidth=0.5, alpha=1, zorder=5)

ax.add_geometries(Reader(fname).geometries(), ccrs.Mercator.GOOGLE, facecolor='none', edgecolor='k', zorder=6, lw=0.4)

#  reduce outline patch linewidths
cbar.outline.set_linewidth(0.4)
ax.outline_patch.set_linewidth(0.4)

# reduce fontsize
cbar.ax.tick_params(labelsize=8) 

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
    
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(8) 

circle = plt.Circle((-70.6828, -33.4450), 0.2, color='k', fill=False, zorder=4, lw=0.5)
ax.add_patch(circle)

plt.savefig('../../../megafires_data/png/CR2MET_1979_2018_trend_full.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()