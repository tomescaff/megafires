import numpy as np
from scipy.optimize import fmin
from scipy.stats import norm
from scipy.stats import genextreme as gev
from scipy.stats import gumbel_r as gum

# maximum likelihood estimation -- norm
def mle_norm_2d(xarr, Tarr, init_params):

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

    xopt = fmin(func=maxlkh, x0 = init_params, args=(xarr, Tarr))
    return xopt

# maximum likelihood estimation -- gev
def mle_gev_2d(xarr, Tarr, init_params):

    def gev_with_trend_shift(x, T, params):
        mu0 = params[0]
        sigma0 = params[1]
        alpha = params[2]
        eta0 = params[3]
        mu = mu0 + alpha*T
        y = gev.pdf(x, eta0, mu, sigma0)
        return y

    def maxlkh(params, *args):
        xs = args[0]
        Ts = args[1]
        f = gev_with_trend_shift
        logp = -sum([np.log(f(x,T, params)) for (x,T) in zip(xs, Ts)])
        return logp

    xopt = fmin(func=maxlkh, x0 = init_params, args=(xarr, Tarr))
    return xopt

# maximum likelihood estimation -- gumbel
def mle_gumbel_2d(xarr, Tarr, init_params):

    def gumbel_with_trend_shift(x, T, params):
        mu0 = params[0]
        sigma0 = params[1]
        alpha = params[2]
        mu = mu0 + alpha*T
        y = gum.pdf(x, mu, sigma0)
        return y

    def maxlkh(params, *args):
        xs = args[0]
        Ts = args[1]
        f = gumbel_with_trend_shift
        logp = -sum([np.log(f(x,T, params)) for (x,T) in zip(xs, Ts)])
        return logp

    xopt = fmin(func=maxlkh, x0 = init_params, args=(xarr, Tarr))
    return xopt

