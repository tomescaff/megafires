import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'

import cartopy.crs as ccrs
import pandas as pd
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

dis_shape = ut.get_regional_shape()

# get CR2MET january tmax series
tmax = ut.get_CR2METv25_jan()
tmax_2017 = tmax.sel(time='2017')

# compute clim
cr2 = ut.get_CR2MET_jan()
cr2 = cr2.sel(time='2017').squeeze()

da = tmax_2017
da = da.squeeze()

mask = ut.get_cl_mask()
da = da*mask

# read stations
stnfile = '../../../megafires_data/series/cr2_tasmax_jan_2017_mean_80.csv'
df = pd.read_csv(stnfile)
stn_lon = df.Longitud.values
stn_lat = df.Latitud.values
stn_val = df.Valor.values

# stnfile = '../../../megafires_data/CR2_explorer/cr2_tasmax_anom2017_refperiod_1991_2020_26deg_40degS.nc'
# ds = xr.open_dataset(stnfile)
# stn_lon = ds.lon.values
# stn_lat = ds.lat.values
# stn_val = ds.data.values

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

ax.add_feature(land, linewidth=0.0, alpha=1, zorder=-2)
ax.add_feature(ocean, alpha = 1, zorder=-1)

# add grid using previous ticks
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='grey', alpha=0.7, linestyle='--', draw_labels=False)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)

# plot the climatology and reshape color bar 
# amwg_blueyellowred
# WhiteBlueGreenYellowRed
cmap=cmaps.WhiteBlueGreenYellowRed
pcm = ax.pcolormesh(da.lon.values, da.lat.values, da.values, cmap=cmap, zorder=1, vmin=15, vmax=35)
cbar = plt.colorbar(pcm, aspect = 40, pad=0.03)

# plot the stations
ax.scatter(stn_lon, stn_lat, s=10, c=stn_val, cmap=cmap, zorder=2, vmin=15, vmax=35, edgecolors='k')

# draw the coastlines
# land = cfeature.NaturalEarthFeature('physical', 'land',  scale=resol, edgecolor='k', facecolor='none')
# ax.add_feature(land, linewidth=0.5, alpha=1, zorder=2)

ax.add_geometries(Reader(fname).geometries(), ccrs.Mercator.GOOGLE, facecolor='none', edgecolor='k', zorder=6, lw=0.4)
# ax.add_geometries([dis_shape], crs=ccrs.PlateCarree(),linewidth=1.5, edgecolor='b',facecolor='none', alpha=1, zorder=7)

#  reduce outline patch linewidths
cbar.outline.set_linewidth(0.4)
ax.spines['geo'].set_linewidth(0.4)

# reduce fontsize

cbar.ax.tick_params(labelsize=8) 

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
    
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(8) 

plt.savefig('../../../megafires_data/png/CR2METv25_jan_2017_actual_value_stations.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()