import sys
import numpy as np

sys.path.append('..')

import processing.gmst as gmst
import processing.lens as lens
import processing.math as pmath
import processing.stations as stns

lens2_gmst_full = gmst.get_gmst_annual_lens2_ensmean()
lens2_tmax_full = lens.get_LENS2_jan_tmax_QNWE()
qn_full = stns.get_QN_tmax_jan()

qn = qn_full.sel(time=slice('1930', '2021'))
lens2_gmst = lens2_gmst_full.sel(time=slice('1930', '2021'))
lens2_tmax = lens2_tmax_full.sel(time=slice('1930', '2021'))

lens2_tmax_anom = lens2_tmax - lens2_tmax.mean('time')
lens2_tmax_corrected = lens2_tmax_anom + qn.mean('time')

lens2_gmst_arr = np.tile(lens2_gmst.values, lens2_tmax_corrected.shape[0])
lens2_tmax_arr = np.ravel(lens2_tmax_corrected.values)

xarr = lens2_tmax_arr
Tarr = lens2_gmst_arr
init_params = [29.57, 1.15, 0.88]

xopt = pmath.mle_norm_2d_fast(xarr, Tarr, init_params)
print(xopt)