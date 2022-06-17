import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import numpy as np
import cmaps

ds = xr.open_dataset('../../../megafires_data/LENS_ALL/LENS_trefht_1979_2006_yearly_trend.nc')
slope = ds['slope'].mean('run')*10

fig = plt.figure(figsize=(12,7))

# define projection
ax = plt.axes(projection=ccrs.PlateCarree())

# set extent of map
# ax.set_extent([-77, -66.5, -56.5, -17], crs=ccrs.PlateCarree())

# define and set  x and y ticks
xticks = np.arange(-180,180+60,60)
yticks = np.arange(-90,90+30,30)
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
pcm = ax.pcolormesh(slope.lon.values, slope.lat.values, slope.values, cmap=cmaps.BlueWhiteOrangeRed, alpha=0.8, zorder=3, vmin=-0.25, vmax=0.25)
cbar = plt.colorbar(pcm, aspect = 60, pad=0.1, location='bottom', orientation='horizontal')

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
cbar.ax.set_title('Mean Temperature Trend 1979-2006 (ºC/dec) ', fontdict={'fontsize':10})

plt.savefig('../../../megafires_data/png/LENS_paper_falvey_garreaud_trefht.png', dpi=300)
plt.show()
