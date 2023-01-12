import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'

import cartopy.crs as ccrs
import pandas as pd
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io.shapereader import Reader
from owslib.wmts import WebMapTileService
import xarray as xr
import numpy as np

import sys
import cmaps

sys.path.append('../../processing')

# da = xr.open_rasterio('../../../megafires_data/geotiff/snapshot-2017-01-20T00_00_00Z.tiff')
da = xr.open_rasterio('../../../megafires_data/geotiff/snapshot-2017-01-26T00_00_00Z.tiff')
# da = xr.open_rasterio('../../../megafires_data/geotiff/snapshot-topo.tiff')
fname = '../../../megafires_data/shp/Regiones/Regional.shp'

# create figure
fig = plt.figure(figsize=(8,7))

# define projection
ax = plt.axes(projection=ccrs.PlateCarree())

# set extent of map
ax.set_extent([-74.5, -68, -40, -26], crs=ccrs.PlateCarree())

# define and set  x and y ticks
xticks = [ -72, -70, -68]
yticks = [ -40, -38, -36, -34, -32, -30, -28, -26]
ax.set_xticks( xticks, crs=ccrs.PlateCarree())
ax.set_yticks( yticks, crs=ccrs.PlateCarree())

# format x and y labels
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

# add grid using previous ticks
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='grey', alpha=0.7, linestyle='--', draw_labels=False)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)


color_tuples = np.array([da.variable.data[0].flatten()/256, da.variable.data[1].flatten()/256, da.variable.data[2].flatten()/256]).transpose()
plt.pcolormesh(da.x, da.y, da.variable.data[0], zorder=1, color=color_tuples)

# draw the coastlines
# land = cfeature.NaturalEarthFeature('physical', 'land',  scale=resol, edgecolor='k', facecolor='none')
# ax.add_feature(land, linewidth=0.5, alpha=1, zorder=2)

ax.add_geometries(Reader(fname).geometries(), ccrs.Mercator.GOOGLE, facecolor='none', edgecolor='k', zorder=6, lw=0.4)
# ax.add_geometries([dis_shape], crs=ccrs.PlateCarree(),linewidth=1.5, edgecolor='b',facecolor='none', alpha=1, zorder=7)

#  reduce outline patch linewidths
ax.spines['geo'].set_linewidth(0.4)

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
    
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(8) 

plt.savefig('../../../megafires_data/png/CR2METv25_satellite_20170126.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()