import sys
import xarray as xr
import numpy as np
from scipy.stats import pearsonr
sys.path.append('../../processing')

import processing.utils as ut
import processing.stations as stns
import processing.fires as fires

qn_jan = stns.get_QN_tmax_jan()
cu_jan = stns.get_CU_tmax_jan()
ch_jan = stns.get_CH_tmax_jan()
qn_feb = stns.get_QN_tmax_feb()
cu_feb = stns.get_CU_tmax_feb()
ch_feb = stns.get_CH_tmax_feb()
regind_jan = stns.get_regindex_tmax_jan()
regind_feb = stns.get_regindex_tmax_feb()

qn_janfeb = stns.get_QN_tmax_janfeb()

####

series = qn_jan.where(qn_jan.time.dt.year != 2017, drop=True)

ff_seas = fires.get_fires_frequency_by_season()
ba_seas = fires.get_burned_area_by_season()
ba_jan = fires.get_burned_area_january()
ba_feb = fires.get_burned_area_february()

ff_seas = ff_seas.where(ff_seas.time.dt.year != 2017, drop=True)
ba_seas = ba_seas.where(ba_seas.time.dt.year != 2017, drop=True)
ba_jan = ba_jan.where(ba_jan.time.dt.year != 2017, drop=True)
ba_feb = ba_feb.where(ba_feb.time.dt.year != 2017, drop=True)


###

series_1964_2021 = series.sel(time=slice('1964', '2021'))
ff_seas_1964_2021 = ff_seas.sel(time=slice('1964', '2021'))
ba_seas_1964_2021 = ba_seas.sel(time=slice('1964', '2021'))

series_1966_2021 = series.sel(time=slice('1966', '2021'))
ff_seas_1966_2021 = ff_seas.sel(time=slice('1966', '2021'))
ba_seas_1966_2021 = ba_seas.sel(time=slice('1966', '2021'))

series_1985_2021 = series.sel(time=slice('1985', '2021'))
ff_seas_1985_2021 = ff_seas.sel(time=slice('1985', '2021'))
ba_seas_1985_2021 = ba_seas.sel(time=slice('1985', '2021'))
ba_jan_1985_2021 = ba_jan.sel(time=slice('1985', '2021'))
ba_feb_1985_2021 = ba_feb.sel(time=slice('1985', '2021'))

series_1995_2021 = series.sel(time=slice('1995', '2021'))
ff_seas_1995_2021 = ff_seas.sel(time=slice('1995', '2021'))
ba_seas_1995_2021 = ba_seas.sel(time=slice('1995', '2021'))
ba_jan_1995_2021 = ba_jan.sel(time=slice('1995', '2021'))
ba_feb_1995_2021 = ba_feb.sel(time=slice('1995', '2021'))

print('1a',pearsonr(series_1964_2021,ff_seas_1964_2021))
print('1b',pearsonr(series_1966_2021,ff_seas_1966_2021))
print('1c',pearsonr(series_1985_2021,ff_seas_1985_2021))
print('1d',pearsonr(series_1995_2021,ff_seas_1995_2021))

print('2a',pearsonr(series_1964_2021,ba_seas_1964_2021))
print('2b',pearsonr(series_1966_2021,ba_seas_1966_2021))
print('2c',pearsonr(series_1985_2021,ba_seas_1985_2021))
print('2d',pearsonr(series_1995_2021,ba_seas_1995_2021))

print('3a',pearsonr(series_1985_2021,ba_jan_1985_2021))
print('3b',pearsonr(series_1995_2021,ba_jan_1995_2021))

print('4a',pearsonr(series_1985_2021,ba_feb_1985_2021))
print('4b',pearsonr(series_1995_2021,ba_feb_1995_2021))
