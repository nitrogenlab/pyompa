from __future__ import division, print_function
import numpy as np


def get_exponential_from_bounds_func(lowerbound, upperbound, alpha, beta):
    assert upperbound > lowerbound, (upperbound, lowerbound)
    return lambda x: beta*(np.exp(
        alpha*np.maximum(0, np.maximum(lowerbound-x, x-upperbound)))-1)


def get_combined_penalty_func(lat_penalty_func, sig0_penalty_func):
    return lambda lat, sig0: lat_penalty_func(lat)+sig0_penalty_func(sig0)


def get_wm_penalty_func(latbounds, sig0bounds):
    LAT_ALPHA = 0.05
    LAT_BETA = 100
    SIG0_ALPHA = 0.5
    SIG0_BETA = 1000
    lat_penalty_func = (get_exponential_from_bounds_func(
        lowerbound=latbounds[0], upperbound=latbounds[1],
        alpha=LAT_ALPHA, beta=LAT_BETA)
        if latbounds is not None else lambda x: 0)
    sig0_penalty_func = (get_exponential_from_bounds_func(
        lowerbound=sig0bounds[0], upperbound=sig0bounds[1],
        alpha=SIG0_ALPHA, beta=SIG0_BETA)
        if sig0bounds is not None else lambda x: 0)
    return get_combined_penalty_func(lat_penalty_func=lat_penalty_func,
                                     sig0_penalty_func=sig0_penalty_func)
