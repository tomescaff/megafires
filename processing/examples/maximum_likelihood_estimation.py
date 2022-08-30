import sys
import numpy as np
from scipy.optimize import fmin
from scipy.stats import norm

sys.path.append('..')

import processing.gmst as gmst
import processing.stations as stns

smt = gmst.get_gmst_annual_5year_smooth().sel(time=slice('1928','2021'))
qn = stns.get_QN_tmax_jan().sel(time=slice('1928','2021'))

def norm_with_trend_shift(x, T, params):

    mu0 = params[0]
    sigma0 = params[1]
    alpha = params[2]

    mu = mu0 + alpha*T
    y = norm.pdf(x, mu, sigma0)
    return y

def maxlkh(params, *args):
    xs = args[0]
    Ts = args[1]
    f = norm_with_trend_shift
    logp = -sum([np.log(f(x,T, params)) for (x,T) in zip(xs, Ts)])
    return logp

tau_i = 15
sigma_i = 2
alpha_i = 0.5

xopt = fmin(func=maxlkh, x0 = [tau_i, sigma_i, alpha_i], args=(qn.values, smt.values))










