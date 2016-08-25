from __future__ import absolute_import

from .stock import Stock

from .markov_multifractal import MarkovSwitchingMultifractal

# We lazy_import the following modules since they import numpy which slows down sage startup
from sage.misc.lazy_import import lazy_import
lazy_import("sage.finance.time_series", ["TimeSeries"])
lazy_import("sage.finance.time_series", ["autoregressive_fit"])
lazy_import("sage.finance.fractal", ["stationary_gaussian_simulation", "fractional_gaussian_noise_simulation", "fractional_brownian_motion_simulation", "multifractal_cascade_random_walk_simulation"])

from .option import black_scholes
