# sage_setup: distribution = sagemath-symbolics
from sage.misc.lazy_import import lazy_import

lazy_import('sage.dynamics.complex_dynamics.mandel_julia',
            ["mandelbrot_plot", "external_ray", "kneading_sequence", "julia_plot"])
del lazy_import
