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

sys.path.append('../../processing')

import processing.utils as ut

# get CR2MET january tmax series
tmax = ut.get_CR2METv25_jan()
tmax_clim = tmax.sel(time=slice('1991', '2020')).mean('time')
tmax_2017 = tmax.sel(time='2017')

da = tmax_2017-tmax_clim
da = da.squeeze()

mask = ut.get_cl_mask()
da = da*mask

# read stations
stnfile = '../../../megafires_data/series/CR2_explorador_tasmax_ene_2017_anom_1991_2020_80_80.csv'
df = pd.read_csv(stnfile)
stn_lon = df.Longitud.values
stn_lat = df.Latitud.values
stn_val = df.Valor.values

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
pcm = ax.pcolormesh(da.lon.values, da.lat.values, da.values, cmap=cmaps.WhiteBlueGreenYellowRed, zorder=1, vmin=0, vmax=3.5)
cbar = plt.colorbar(pcm, aspect = 40, pad=0.03)

# plot the stations
ax.scatter(stn_lon, stn_lat, s=10, c=stn_val, cmap=cmaps.WhiteBlueGreenYellowRed, zorder=2, vmin=0, vmax=3.5, edgecolors='k')

# draw the coastlines
#land = cfeature.NaturalEarthFeature('physical', 'land',  scale=resol, edgecolor='k', facecolor='none')
#ax.add_feature(land, linewidth=0.5, alpha=1, zorder=2)

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

# circle_qn = plt.Circle((-70.6828, -33.4450), 0.2, color='k', fill=False, zorder=4, lw=0.5)
# circle_cu = plt.Circle((-71.2167, -34.9664), 0.2, color='k', fill=False, zorder=4, lw=0.5)
# circle_ch = plt.Circle((-72.0400, -36.5872), 0.2, color='k', fill=False, zorder=4, lw=0.5)
# circle_cc = plt.Circle((-73.0622, -36.7792), 0.2, color='k', fill=False, zorder=4, lw=0.5)

# ax.add_patch(circle_qn)
# ax.add_patch(circle_cu)
# ax.add_patch(circle_ch)
# ax.add_patch(circle_cc)

plt.savefig('../../../megafires_data/png/CR2METv25_jan_2017_tmax_anom_full_stations.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()