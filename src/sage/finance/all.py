from sage.misc.lazy_import import lazy_import

lazy_import('sage.finance.stock', 'Stock')
lazy_import('sage.finance.markov_multifractal', 'MarkovSwitchingMultifractal')

# We lazy_import the following modules since they import numpy which
# slows downsage startup
lazy_import('sage.finance.time_series', 'TimeSeries')
lazy_import('sage.finance.time_series', 'autoregressive_fit')
lazy_import('sage.finance.fractal',
            ['stationary_gaussian_simulation',
             'fractional_gaussian_noise_simulation',
             'fractional_brownian_motion_simulation',
             'multifractal_cascade_random_walk_simulation'])
lazy_import('sage.finance.option', 'black_scholes')
