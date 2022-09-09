import sys
import numpy as np

sys.path.append('..')

import processing.gmst as gmst
import processing.lens as lens
import processing.math as pmath
import processing.stations as stns

lens1_gmst_full = gmst.get_gmst_annual_lens1_ensmean()
lens1_tmax_full = lens.get_LENS_jan_tmax_QNWE()
qn_full = stns.get_QN_tmax_jan()

qn = qn_full.sel(time=slice('1930', '2021'))
lens1_gmst = lens1_gmst_full.sel(time=slice('1930', '2021'))
lens1_tmax = lens1_tmax_full.sel(time=slice('1930', '2021'))

lens1_tmax_anom = lens1_tmax - lens1_tmax.mean('time')
lens1_tmax_corrected = lens1_tmax_anom + qn.mean('time')

lens1_gmst_arr = np.tile(lens1_gmst.values, lens1_tmax_corrected.shape[0])
lens1_tmax_arr = np.ravel(lens1_tmax_corrected.values)

xarr = lens1_tmax_arr
Tarr = lens1_gmst_arr
init_params = [29.54, 1.03, 1.11]

xopt = pmath.mle_norm_2d_fast(xarr, Tarr, init_params)
print(xopt)