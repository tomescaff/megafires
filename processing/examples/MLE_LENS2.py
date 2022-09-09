import sys
import numpy as np

sys.path.append('..')

import processing.gmst as gmst
import processing.lens as lens
import processing.math as pmath

lens2_gmst = gmst.get_gmst_annual_lens2_ensmean()
lens2_tmax = lens.get_LENS2_jan_tmax_QNWE()

lens2_gmst_arr = np.tile(lens2_gmst.values, lens2_tmax.shape[0])
lens2_tmax_arr = np.ravel(lens2_tmax.values)

xarr = lens2_tmax_arr
Tarr = lens2_gmst_arr
init_params = [15, 2, 0.5]

xopt = pmath.mle_norm_2d_fast(xarr, Tarr, init_params)