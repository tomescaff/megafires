import matplotlib
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
import xarray as xr

sys.path.append('../../processing')

import processing.utils as ut

# get CR2MET january tmax return period

da = xr.open_dataset('../../../megafires_data/output/ee2017_return_period.nc')['return_period']

lats = [-33.45, -34.97, -36.59, -29.92, -36.78, -38.77, -39.65, -33.39, -33.45, -33.66]
lons = [-70.68, -71.22, -72.04, -71.20, -73.06, -72.64, -73.08, -70.79, -70.55, -71.61]
retp = [1.4e3, 2.4e2, 1.0e1, 2.6e1, 5.6e1, 1.6e0, 1.6e0, 2.3e3, 1.3e3, 2.7e1]

fname = '../../../megafires_data/shp/Regiones/Regional.shp'

# create figure
fig = plt.figure(figsize=(8,7))

# define projection
ax = plt.axes(projection=ccrs.PlateCarree())

# set extent of map
ax.set_extent([-74.5, -68, -40, -26], crs=ccrs.PlateCarree())

# define and set  x and y ticks
xticks = [ -74, -72, -70, -68]
yticks = [ -40, -38, -36, -34, -32, -30, -28, -26]
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
pcm = ax.pcolormesh(da.lon.values, da.lat.values, da.values, cmap=cmaps.WhiteBlueGreenYellowRed, zorder=4, vmin=1, vmax=1e6, norm=matplotlib.colors.LogNorm())
cbar = plt.colorbar(pcm, aspect = 40, pad=0.03)

# plot the stations
ax.scatter(lons, lats, s=10, c=retp, cmap=cmaps.WhiteBlueGreenYellowRed, zorder=5, edgecolors='k', norm=matplotlib.colors.LogNorm(vmin=1, vmax=1e6))

# draw the coastlines
# land = cfeature.NaturalEarthFeature('physical', 'land',  scale=resol, edgecolor='k', facecolor='none')
# ax.add_feature(land, linewidth=0.5, alpha=1, zorder=5)

ax.add_geometries(Reader(fname).geometries(), ccrs.Mercator.GOOGLE, facecolor='none', edgecolor='k', zorder=6, lw=0.4)

#  reduce outline patch linewidths
cbar.outline.set_linewidth(0.4)
ax.spines['geo'].set_linewidth(0.4)

# reduce fontsize
cbar.ax.tick_params(labelsize=8) 

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
    
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(8) 

#circle_qn = plt.Circle((-70.6828, -33.4450), 0.2, color='k', fill=False, zorder=7, lw=0.5)
#circle_cu = plt.Circle((-71.2167, -34.9664), 0.2, color='k', fill=False, zorder=7, lw=0.5)
#circle_ch = plt.Circle((-72.0400, -36.5872), 0.2, color='k', fill=False, zorder=7, lw=0.5)
#circle_cc = plt.Circle((-73.0622, -36.7792), 0.2, color='k', fill=False, zorder=7, lw=0.5)

plt.savefig('../../../megafires_data/png/CR2METv25_jan_1960_2021_return_period_2017_stations.png', dpi=300)#, bbox_inches = 'tight', pad_inches = 0)
plt.show()