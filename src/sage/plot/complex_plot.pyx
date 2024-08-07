"""
Complex plots

AUTHORS:

- Robert Bradshaw (2009): initial version
- David Lowry-Duda (2022): incorporate matplotlib colormaps
"""
# ****************************************************************************
#       Copyright (C) 2009 Robert Bradshaw <robertwb@math.washington.edu>,
#       Copyright (C) 2022 David Lowry-Duda <david@lowryduda.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off, sig_check


# Note: don't import matplotlib at module level. It takes a surprisingly
#       long time to load!
cimport numpy as cnumpy

from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options
from sage.rings.complex_double cimport ComplexDoubleElement
from sage.arith.srange import srange

from libc.math cimport hypot, atan2, atan, log, pow, sqrt
from sage.arith.constants cimport M_PI as PI

from sage.libs.gsl.complex cimport *


DEFAULT_LOGARITHMIC_CONTOUR_BASE = 2
DEFAULT_LINEAR_CONTOUR_BASE = 10


cdef inline ComplexDoubleElement new_CDF_element(double x, double y):
    z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
    GSL_SET_COMPLEX(&z._complex, x, y)
    return z


cdef inline double mag_to_lightness(double r, double rate=0.5) noexcept:
    """
    Return a lightness for the given magnitude.

    Small magnitudes are darker and large magnitudes are lighter, and the
    lightness smoothly transitions.

    Tweak the rate to adjust how the magnitude affects the color. For
    instance, changing ``rate`` to `1` will cause anything near zero to be
    much darker and poles to be much lighter, while `0.25` would cause the
    reverse effect.

    INPUT:

    - ``r`` -- a nonnegative real number; the magnitude

    - ``rate`` -- a positive real number; how quickly changes in magnitude
      affect changes in output lightness

    OUTPUT:

    A value between `-1` (black) and `+1` (white), inclusive.

    .. SEEALSO::

        - :func:`sage.plot.complex_plot.cyclic_logarithmic_mag_to_lightness`
        - :func:`sage.plot.complex_plot.mag_and_arg_to_lightness`

    EXAMPLES:

    This tests it implicitly::

        sage: from sage.plot.complex_plot import complex_to_rgb
        sage: complex_to_rgb([[0, 1, 10]])  # abs tol 1e-4
        array([[[0.        , 0.        , 0.        ],
                [0.77172568, 0.        , 0.        ],
                [1.        , 0.22134776, 0.22134776]]])

    Changing ``rate`` changes the rate of growth::

        sage: complex_to_rgb([[0, 1, 10]], dark_rate=1)  # abs tol 1e-4
        array([[[0.        , 0.        , 0.        ],
                [0.77172568, 0.        , 0.        ],
                [1.        , 0.49693961, 0.49693961]]])
    """
    if rate == 0.5:
        return atan(log(sqrt(r)+1)) * (4/PI) - 1
    else:
        return atan(log(pow(r, rate)+1)) * (4/PI) - 1


cdef inline double cyclic_logarithmic_mag_to_lightness(double r, double base=2) noexcept:
    r"""
    Return a lightness for the given magnitude.

    This modifies the lightness around magnitudes of size `base^n` for integer
    `n` to create the appearance of logarithmically spaced contours.

    INPUT:

    - ``r`` -- a nonnegative real number; the magnitude

    - ``base`` -- a real number (default: `2`); contours will appear at integer
      powers of ``base``. This should be greater than `1`.

    OUTPUT:

    A value between `-1` (black) and `+1` (white), inclusive.

    .. SEEALSO::

        - :func:`sage.plot.complex_plot.mag_to_lightness`
        - :func:`sage.plot.complex_plot.mag_and_arg_to_lightness`
        - :func:`sage.plot.complex_plot.cyclic_linear_mag_to_lightness`

    EXAMPLES:

    This tests it implicitly::

        sage: from sage.plot.complex_plot import complex_to_rgb
        sage: complex_to_rgb([[0, 1, 10]], contoured=True,  # abs tol 1e-4
        ....:                contour_type='logarithmic')
        array([[[1.        , 0.        , 0.        ],
                [1.        , 0.15      , 0.15      ],
                [0.98903595, 0.        , 0.        ]]])

    We set contours to be multiples of `5` apart::

        sage: complex_to_rgb([[0, 1, 10]], contoured=True,  # abs tol 1e-4
        ....:                contour_type='logarithmic', contour_base=5)
        array([[[1.        , 0.        , 0.        ],
                [1.        , 0.15      , 0.15      ],
                [0.93466172, 0.        , 0.        ]]])
    """
    if r < 1e-10:
        return 0.0
    rem  = (log(r) / log(base)) % 1
    if rem < 0:  # Choose positive mod representative
        rem += 1
    return .15 - rem/2.


cdef inline double cyclic_linear_mag_to_lightness(double r, double base=10) noexcept:
    r"""
    Return a lightness for the given magnitude.

    This modifies the lightness around magnitudes of size `n*base` for integer
    `n` to create the appearance of linearly spaced contours.

    INPUT:

    - ``r`` -- a nonnegative real number; the magnitude

    - ``base`` -- a positive real number (default: `10`); contours will appear
      at integer multiples of ``base``

    OUTPUT:

    A value between `-1` (black) and `+1` (white), inclusive.

    .. SEEALSO::

        - :func:`sage.plot.complex_plot.mag_to_lightness`
        - :func:`sage.plot.complex_plot.mag_and_arg_to_lightness`
        - :func:`sage.plot.complex_plot.cyclic_logarithmic_mag_to_lightness`

    EXAMPLES:

    This tests it implicitly::

        sage: from sage.plot.complex_plot import complex_to_rgb
        sage: complex_to_rgb([[1, 5, 11]], contoured=True,  # abs tol 1e-4
        ....:                contour_type='linear')
        array([[[1. , 0.1, 0.1],
                [0.9, 0. , 0. ],
                [1. , 0.1, 0.1]]])

    In the above example, note that both `1` and `11` have the same imaginary
    part and are precisely `10` (the default contour separation) apart. If we
    set contours to be multiples of `3` apart, the values are no longer the
    same, but the values for `5` and `11` should be::

        sage: complex_to_rgb([[1, 5, 11]], contoured=True,  # abs tol 1e-4
        ....:                contour_type='linear', contour_base=3)
        array([[[0.98333333, 0.        , 0.        ],
                [0.81666667, 0.        , 0.        ],
                [0.81666667, 0.        , 0.        ]]])
    """
    rem = (r / base) % 1
    if rem < 0:  # Choose positive mod representative
        rem += 1
    return .15 - rem/2.


cdef inline double mag_and_arg_to_lightness(double r, double arg,
                                            double base=2, int nphases=10) noexcept:
    r"""
    Return a lightness for the given magnitude and argument.

    This modifies the lightness around magnitudes of size ``base^n`` for
    integer `n`, and also around arguments of the form `(2 \pi) / nphases * n` for
    integer `n`. This creates the appearance of tiles, bounded by
    (base ``base``) logarithmically spaced magnitude contours and by
    (``nphases`` many) evenly spaced argument contours.

    INPUT:

    - ``r`` -- a nonnegative real number

    - ``arg`` -- a real number

    - ``base`` -- a positive real number (default: ``2``); contours will appear
      at integer powers of ``base``. This should be greater than `1`.

    - ``nphases`` -- positive integer (default: `10`); how many phase
      contours to represent a change of argument by `2 \pi`

    OUTPUT:

    A value between `-1` (black) and `+1` (white), inclusive.

    .. SEEALSO::

        - :func:`sage.plot.complex_plot.mag_to_lightness`
        - :func:`sage.plot.complex_plot.cyclic_logarithmic_mag_to_lightness`

    EXAMPLES:

    This tests it implicitly::

        sage: from sage.plot.complex_plot import complex_to_rgb
        sage: complex_to_rgb([[0, 1, 10]], tiled=True)  # abs tol 1e-4
        array([[[1.        , 0.        , 0.        ],
                [1.        , 0.15      , 0.15      ],
                [1.        , 0.06951798, 0.06951798]]])
        sage: complex_to_rgb([[0, 1 + 1j, -3 - 5j]], tiled=True)  # abs tol 1e-4
        array([[[1.        , 0.        , 0.        ],
                [0.9625    , 0.721875  , 0.        ],
                [0.        , 0.01371897, 0.85409323]]])

    Adjusting the tiling parameters should create relatively small
    differences::

        sage: complex_to_rgb([[0, 1 + 1j, -3 - 5j]],  # abs tol 1e-4
        ....:                tiled=True, nphases=15)
        array([[[1.        , 0.        , 0.        ],
                [0.80625   , 0.6046875 , 0.        ],
                [0.        , 0.01243417, 0.77410628]]])
        sage: complex_to_rgb([[0, 1 + 1j, -3 - 5j]],  # abs tol 1e-4
        ....:                tiled=True, contour_base=5, nphases=15)
        array([[[1.        , 0.        , 0.        ],
                [0.87741543, 0.65806157, 0.        ],
                [0.        , 0.01423401, 0.88615776]]])
    """
    if r < 1e-10:
        return 0.0
    cdef double r_rem, arg_rem
    r_rem  = (log(r) / log(base)) % 1
    arg_rem = (nphases * arg / (2*PI)) % 1
    if r_rem < 0:  # Choose positive mod representatives
        r_rem += 1
    if arg_rem < 0:
        arg_rem += 1
    return 0.15 - (r_rem)/4. - (arg_rem)/4.


def complex_to_rgb(z_values, contoured=False, tiled=False,
                   contour_type='logarithmic', contour_base=None,
                   dark_rate=0.5, nphases=10):
    r"""
    Convert a grid of complex numbers to a grid of rgb values using a default
    choice of colors.

    INPUT:

    - ``z_values`` -- a grid of complex numbers, as a list of lists

    - ``contoured`` -- boolean (default: ``False``); causes magnitude to be
      indicated through contour-like adjustments to lightness

    - ``tiled`` -- boolean (default: ``False``); causes magnitude and argument to
      be indicated through contour-like adjustments to lightness

    - ``nphases`` -- positive integer (default: 10); when ``tiled=True``,
      this is the number of divisions the phase is divided into

    - ``contour_type`` -- either ``'logarithmic'``, or ``'linear'`` (default:
      ``'logarithmic'``); causes added contours to be of given type when
      ``contoured=True``.

    - ``contour_base`` -- positive integer; when ``contour_type`` is
      ``'logarithmic'``, this sets logarithmic contours at multiples of
      ``contour_base`` apart. When ``contour_type`` is ``'linear'``, this sets
      contours at distances of ``contour_base`` apart. If ``None``, then a
      default is chosen depending on ``contour_type``.

    - ``dark_rate`` -- a positive number (default: 0.5); affects how quickly
      magnitudes affect how light/dark the image is. When there are contours,
      this affects how visible each contour is. Large values (near `1.0`) have
      very strong, immediate effects, while small values (near `0.0`) have
      gradual effects.

    OUTPUT:

    An `N \times M \times 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r,g,b) tuple.

    .. SEEALSO::

        :func:`sage.plot.complex_plot.complex_to_cmap_rgb`

    EXAMPLES:

    We can call this on grids of complex numbers::

        sage: from sage.plot.complex_plot import complex_to_rgb
        sage: complex_to_rgb([[0, 1, 1000]])  # abs tol 1e-4
        array([[[0.        , 0.        , 0.        ],
                [0.77172568, 0.        , 0.        ],
                [1.        , 0.64421177, 0.64421177]]])
        sage: complex_to_rgb([[0, 1j, 1000j]])  # abs tol 1e-4
        array([[[0.        , 0.        , 0.        ],
                [0.38586284, 0.77172568, 0.        ],
                [0.82210588, 1.        , 0.64421177]]])
        sage: complex_to_rgb([[0, 1, 1000]], contoured=True)   # abs tol 1e-4
        array([[[1.        , 0.        , 0.        ],
                [1.        , 0.15      , 0.15      ],
                [0.66710786, 0.        , 0.        ]]])
        sage: complex_to_rgb([[0, 1, 1000]], tiled=True)   # abs tol 1e-4
        array([[[1.        , 0.        , 0.        ],
                [1.        , 0.15      , 0.15      ],
                [0.90855393, 0.        , 0.        ]]])

    We can change contour types and the distances between contours::

        sage: complex_to_rgb([[0, 1 + 1j, 3 + 4j]],  # abs tol 1e-4
        ....:                contoured=True, contour_type='logarithmic', contour_base=3)
        array([[[1.        , 0.        , 0.        ],
                [0.99226756, 0.74420067, 0.        ],
                [0.91751324, 0.81245954, 0.        ]]])
        sage: complex_to_rgb([[0, 1 + 1j, 3 + 4j]],  # abs tol 1e-4
        ....:                contoured=True, contour_type='linear', contour_base=3)
        array([[[1.        , 0.15      , 0.15      ],
                [0.91429774, 0.6857233 , 0.        ],
                [0.81666667, 0.72315973, 0.        ]]])

    Lowering ``dark_rate`` causes colors to go to black more slowly near `0`::

        sage: complex_to_rgb([[0, 0.5, 1]], dark_rate=0.4)  # abs tol 1e-4
        array([[[0.        , 0.        , 0.        ],
                [0.65393731, 0.        , 0.        ],
                [0.77172568, 0.        , 0.        ]]])
        sage: complex_to_rgb([[0, 0.5, 1]], dark_rate=0.2)  # abs tol 1e-4
        array([[[0.        , 0.        , 0.        ],
                [0.71235886, 0.        , 0.        ],
                [0.77172568, 0.        , 0.        ]]])
    """
    import numpy as np

    cdef unsigned int i, j, imax, jmax
    cdef double x, y, mag, arg
    cdef double lightness, hue, top, bot
    cdef double r, g, b
    cdef int ihue
    cdef ComplexDoubleElement z
    from sage.rings.complex_double import CDF

    imax = len(z_values)
    jmax = len(z_values[0])
    cdef cnumpy.ndarray[cnumpy.float_t, ndim=3, mode='c'] rgb = np.empty(dtype=float, shape=(imax, jmax, 3))

    if contour_base is None:
        if contour_type == "linear":
            contour_base = DEFAULT_LINEAR_CONTOUR_BASE
        else:
            contour_base = DEFAULT_LOGARITHMIC_CONTOUR_BASE
    if contour_type not in ("linear", "logarithmic"):
        raise ValueError("Unrecognized contour type argument {}.".format(contour_type))
    if contour_base <= 0:
        raise ValueError("contour_base must be positive")

    sig_on()
    for i in range(imax):
        row = z_values[i]
        for j in range(jmax):
            zz = row[j]
            if type(zz) is ComplexDoubleElement:
                z = <ComplexDoubleElement>zz
            else:
                z = CDF(zz)
            x = GSL_REAL(z._complex)
            y = GSL_IMAG(z._complex)
            mag = hypot(x, y)
            arg = atan2(y, x) # math module arctan has range from -pi to pi, so cut along negative x-axis

            if tiled:
                lightness = mag_and_arg_to_lightness(
                    mag, arg, base=contour_base, nphases=nphases
                )
            elif contoured:
                if contour_type == "logarithmic":
                    lightness = cyclic_logarithmic_mag_to_lightness(mag, base=contour_base)
                else:
                    lightness = cyclic_linear_mag_to_lightness(mag, base=contour_base)
            else:
                lightness = mag_to_lightness(mag, rate=dark_rate)

            if lightness < 0: # in hsv, variable value, full saturation (s=1, v=1+lightness)
                bot = 0
                top = (1+lightness)
            else: # in hsv, variable saturation, full value (v=1, s=1-lightness)
                bot = lightness
                top = 1

            # Note that does same thing as colorsys module hsv_to_rgb
            # for this setup, but in Cython
            hue = 3*arg/PI

            if hue < 0:
                hue += 6  # usual hsv hue is thus h=arg/(2*pi) for positive, h=arg/(2*PI)+1 for negative
            ihue = <int>hue
            if ihue == 0:
                r = top
                g = bot + hue * (top-bot)
                b = bot
            elif ihue == 1:
                r = bot + (2-hue) * (top-bot)
                g = top
                b = bot
            elif ihue == 2:
                r = bot
                g = top
                b = bot + (hue-2) * (top-bot)
            elif ihue == 3:
                r = bot
                g = bot + (4-hue) * (top-bot)
                b = top
            elif ihue == 4:
                r = bot + (hue-4) * (top-bot)
                g = bot
                b = top
            else:
                r = top
                g = bot
                b = bot + (6-hue) * (top-bot)

            rgb[i, j, 0] = r
            rgb[i, j, 1] = g
            rgb[i, j, 2] = b

    sig_off()
    nan_indices = np.isnan(rgb).any(-1)     # Mask for undefined points
    rgb[nan_indices] = 1                    # Make nan_indices white
    return rgb


def complex_to_cmap_rgb(z_values, cmap='turbo', contoured=False, tiled=False,
                        contour_type='logarithmic', contour_base=None,
                        dark_rate=0.5, nphases=10):
    r"""
    Convert a grid of complex numbers to a grid of rgb values using colors
    taken from given colormap.

    INPUT:

    - ``z_values`` -- a grid of complex numbers, as a list of lists

    - ``cmap`` -- the string name of a matplotlib colormap, or an instance
      of a matplotlib Colormap (default: ``'turbo'``)

    - ``contoured`` -- boolean (default: ``False``); causes magnitude to be
      indicated through contour-like adjustments to lightness

    - ``tiled`` -- boolean (default: ``False``); causes magnitude and argument to
      be indicated through contour-like adjustments to lightness

    - ``nphases`` -- positive integer (default: 10); when ``tiled=True``,
      this is the number of divisions the phase is divided into

    - ``contour_type`` -- either ``'logarithmic'``, or ``'linear'`` (default:
      ``'logarithmic'``); causes added contours to be of given type when
      ``contoured=True``.

    - ``contour_base`` -- positive integer; when ``contour_type`` is
      ``'logarithmic'``, this sets logarithmic contours at multiples of
      ``contour_base`` apart. When ``contour_type`` is ``'linear'``, this sets
      contours at distances of ``contour_base`` apart. If ``None``, then a
      default is chosen depending on ``contour_type``.

    - ``dark_rate`` -- a positive number (default: 0.5); affects how quickly
      magnitudes affect how light/dark the image is. When there are contours,
      this affects how visible each contour is. Large values (near 1.0) have
      very strong, immediate effects, while small values (near 0.0) have
      gradual effects.

    OUTPUT:

    An `N \times M \times 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r, g, b) tuple.

    .. SEEALSO::

        :func:`sage.plot.complex_plot.complex_to_rgb`

    EXAMPLES:

    We can call this on grids of complex numbers::

        sage: from sage.plot.complex_plot import complex_to_cmap_rgb
        sage: complex_to_cmap_rgb([[0, 1, 1000]])  # abs tol 1e-4
        array([[[0.        , 0.        , 0.        ],
                [0.49669808, 0.76400071, 0.18024425],
                [0.87320419, 0.99643856, 0.72730967]]])
        sage: complex_to_cmap_rgb([[0, 1, 1000]], cmap='viridis')  # abs tol 1e-4
        array([[[0.        , 0.        , 0.        ],
                [0.0984475 , 0.4375291 , 0.42487821],
                [0.68959896, 0.84592555, 0.84009311]]])


    We can change contour types and the distances between contours::

        sage: complex_to_cmap_rgb([[0, 1 + 1j, 3 + 4j]], contoured=True,  # abs tol 1e-4
        ....:                     contour_type='logarithmic', contour_base=3)
        array([[[0.64362   , 0.98999   , 0.23356   ],
                [0.93239357, 0.81063338, 0.21955399],
                [0.95647342, 0.74861225, 0.14963982]]])
        sage: complex_to_cmap_rgb([[0, 1 + 1j, 3 + 4j]], cmap='turbo',   # abs tol 1e-4
        ....:                     contoured=True, contour_type='linear', contour_base=3)
        array([[[0.71246796, 0.9919238 , 0.3816262 ],
                [0.92617785, 0.79322304, 0.14779989],
                [0.95156284, 0.72025117, 0.05370383]]])

    We see that changing ``dark_rate`` affects how visible contours are. In this
    example, we set ``contour_base=5`` and note that the points `0` and `1 + i`
    are far away from contours, but `2.9 + 4i` is near (and just below) a
    contour. Raising ``dark_rate`` should have strong effects on the last
    coloration and weaker effects on the others::

        sage: complex_to_cmap_rgb([[0, 1 + 1j, 2.9 + 4j]], cmap='turbo',  # abs tol 1e-4
        ....:                     contoured=True, dark_rate=0.05, contour_base=5)
        array([[[0.64362   , 0.98999   , 0.23356   ],
                [0.93334746, 0.81330523, 0.23056563],
                [0.96357185, 0.75337736, 0.19440913]]])
        sage: complex_to_cmap_rgb([[0, 1 + 1j, 2.9 + 4j]], cmap='turbo',  # abs tol 1e-4
        ....:                     contoured=True, dark_rate=0.85, contour_base=5)
        array([[[0.64362   , 0.98999   , 0.23356   ],
                [0.93874682, 0.82842892, 0.29289564],
                [0.57778954, 0.42703289, 0.02612716]]])
    """
    import numpy as np
    import matplotlib as mpl

    if isinstance(cmap, str):
        cmap = mpl.colormaps[cmap]

    if contour_base is None:
        if contour_type == "linear":
            contour_base = DEFAULT_LINEAR_CONTOUR_BASE
        else:
            contour_base = DEFAULT_LOGARITHMIC_CONTOUR_BASE
    if contour_type not in ("linear", "logarithmic"):
        raise ValueError("Unrecognized contour type argument {}.".format(contour_type))
    if contour_base <= 0:
        raise ValueError("contour_base must be positive")

    cdef unsigned int i, j, imax, jmax
    cdef double x, y, mag, arg
    cdef double lightness_delta
    cdef ComplexDoubleElement z
    from sage.rings.complex_double import CDF

    imax = len(z_values)
    jmax = len(z_values[0])
    cdef cnumpy.ndarray[cnumpy.float_t, ndim=3, mode='c'] als = np.empty(dtype=float, shape=(imax, jmax, 2))

    sig_on()
    for i in range(imax):
        row = z_values[i]
        for j in range(jmax):
            zz = row[j]
            if type(zz) is ComplexDoubleElement:
                z = <ComplexDoubleElement>zz
            else:
                z = CDF(zz)
            x = GSL_REAL(z._complex)
            y = GSL_IMAG(z._complex)
            mag = hypot(x, y)
            arg = atan2(y, x) # math module arctan has range from -pi to pi, so cut along negative x-axis
            if tiled:
                lightness_delta = mag_and_arg_to_lightness(
                    mag, arg, base=contour_base, nphases=nphases
                )
            elif contoured:
                if contour_type == "logarithmic":
                    lightness_delta = cyclic_logarithmic_mag_to_lightness(mag, base=contour_base)
                else:
                    lightness_delta = cyclic_linear_mag_to_lightness(mag, base=contour_base)
            else:
                lightness_delta = mag_to_lightness(mag, rate=dark_rate)
            als[i, j, 0] = arg
            als[i, j, 1] = lightness_delta
    sig_off()

    args = als[:,:,0]
    nan_indices = np.isnan(als).any(-1)            # Mask for undefined points
    normalized_colors = cmap((args + PI) / (2 * PI)) # break on negative reals
    normalized_colors = normalized_colors[:,:,:3]  # discard alpha channel
    lightdeltas = als[:,:,1]

    if tiled or contoured:
        rgbs = add_contours_to_rgb(normalized_colors, lightdeltas, dark_rate=dark_rate)
    else:
        rgbs = add_lightness_smoothing_to_rgb(normalized_colors, lightdeltas)

    # Apply mask, making nan_indices white
    rgbs[nan_indices] = 1

    return rgbs


def add_lightness_smoothing_to_rgb(rgb, delta):
    r"""
    Return an rgb array from given array of colors and lightness adjustments.

    This smoothly adds lightness from black (when ``delta`` is `-1`) to white
    (when ``delta`` is `1`).

    Each input `(r, g, b)` is modified by ``delta`` to be lighter or darker
    depending on the size of ``delta``. When ``delta`` is `-1`, the output is
    black. When ``delta`` is `+1`, the output is white. Colors
    piecewise-linearly vary from black to the initial `(r, g, b)` to white.

    We assume that the ``delta`` values come from a function like
    :func:`sage.plot.complex_plot.mag_to_lightness`, which maps magnitudes to
    the range `[-1, +1]`.

    INPUT:

    - ``rgb`` -- a grid of length 3 tuples `(r, g, b)`, as an
      `N \times M \times 3` numpy array

    - ``delta`` -- a grid of values as an `N \times M` numpy array; these
      represent how much to change the lightness of each `(r, g, b)`. Values
      should be in `[-1, 1]`.

    OUTPUT:

    An `N \times M \times 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r, g, b) tuple.

    .. SEEALSO::

        - :func:`sage.plot.complex_plot.complex_to_rgb`
        - :func:`sage.plot.complex_plot.add_contours_to_rgb`

    EXAMPLES:

    We can call this on grids of values::

        sage: # needs numpy
        sage: import numpy as np
        sage: from sage.plot.complex_plot import add_lightness_smoothing_to_rgb
        sage: add_lightness_smoothing_to_rgb(  # abs tol 1e-4
        ....:     np.array([[[0, 0.25, 0.5]]]), np.array([[0.75]]))
        array([[[0.75  , 0.8125, 0.875 ]]])
        sage: add_lightness_smoothing_to_rgb(  # abs tol 1e-4
        ....:     np.array([[[0, 0.25, 0.5]]]), np.array([[0.75]]))
        array([[[0.75  , 0.8125, 0.875 ]]])
    """
    import numpy as np
    delta = delta[:,:,np.newaxis]
    delta_pos = delta > 0.0
    rgb = (1.0 - np.abs(delta))*(rgb - delta_pos) + delta_pos
    rgb = np.clip(rgb, 0.0, 1.0)
    return rgb


def add_contours_to_rgb(rgb, delta, dark_rate=0.5):
    r"""
    Return an rgb array from given array of `(r, g, b)` and `(delta)`.

    Each input `(r, g, b)` is modified by ``delta`` to be lighter or darker
    depending on the size of ``delta``. Negative ``delta`` values darken the
    color, while positive ``delta`` values lighten the pixel.

    We assume that the ``delta`` values come from a function like
    :func:`sage.plot.complex_plot.mag_to_lightness`, which maps magnitudes to
    the range `[-1, +1]`.

    INPUT:

    - ``rgb`` -- a grid of length 3 tuples `(r, g, b)`, as an `N \times M
      \times 3` numpy array

    - ``delta`` -- a grid of values as an `N \times M` numpy array; these
      represent how much to change the lightness of each `(r, g, b)`. Values
      should be in `[-1, 1]`.

    - ``dark_rate`` -- a positive number (default: 0.5); affects how
      strongly visible the contours appear

    OUTPUT:

    An `N \times M \times 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r, g, b) tuple.

    .. SEEALSO::

        - :func:`sage.plot.complex_plot.complex_to_rgb`,
        - :func:`sage.plot.complex_plot.add_lightness_smoothing_to_rgb`

    ALGORITHM:

    Each pixel and lightness-delta is mapped from `(r, g, b, delta) \mapsto
    (h, l, s, delta)` using the standard RGB-to-HLS formula.

    Then the lightness is adjusted via `l \mapsto l' = l + 0.5 \cdot delta`.

    Finally map `(h, l', s) \mapsto (r, g, b)` using the standard HLS-to-RGB
    formula.

    EXAMPLES::

        sage: # needs numpy
        sage: import numpy as np
        sage: from sage.plot.complex_plot import add_contours_to_rgb
        sage: add_contours_to_rgb(np.array([[[0, 0.25, 0.5]]]),  # abs tol 1e-4
        ....:                     np.array([[0.75]]))
        array([[[0.25 , 0.625, 1.   ]]])
        sage: add_contours_to_rgb(np.array([[[0, 0, 0]]]),  # abs tol 1e-4
        ....:                     np.array([[1]]))
        array([[[0.5, 0.5, 0.5]]])
        sage: add_contours_to_rgb(np.array([[[1, 1, 1]]]),  # abs tol 1e-4
        ....:                     np.array([[-0.5]]))
        array([[[0.75, 0.75, 0.75]]])

    Raising ``dark_rate`` leads to bigger adjustments::

        sage: add_contours_to_rgb(np.array([[[0.5, 0.5, 0.5]]]),  # abs tol 1e-4        # needs numpy
        ....:                     np.array([[0.5]]), dark_rate=0.1)
        array([[[0.55, 0.55, 0.55]]])
        sage: add_contours_to_rgb(np.array([[[0.5, 0.5, 0.5]]]),  # abs tol 1e-4        # needs numpy
        ....:                     np.array([[0.5]]), dark_rate=0.5)
        array([[[0.75, 0.75, 0.75]]])
    """
    import numpy as np
    hls = rgb_to_hls(rgb)
    hls[..., 1] += dark_rate * delta
    hls = np.clip(hls, 0.0, 1.0)
    rgb = hls_to_rgb(hls)
    return rgb


class ComplexPlot(GraphicPrimitive):
    """
    The GraphicsPrimitive to display complex functions in using the domain
    coloring method

    INPUT:

    - ``rgb_data`` -- an array of colored points to be plotted

    - ``x_range`` -- a minimum and maximum x value for the plot

    - ``y_range`` -- a minimum and maximum y value for the plot

    TESTS::

        sage: p = complex_plot(lambda z: z^2-1, (-2, 2), (-2, 2))
    """
    def __init__(self, rgb_data, x_range, y_range, options):
        """
        TESTS::

            sage: p = complex_plot(lambda z: z^2-1, (-2, 2), (-2, 2))
        """
        self.x_range = x_range
        self.y_range = y_range
        self.x_count = len(rgb_data)
        self.y_count = len(rgb_data[0])
        self.rgb_data = rgb_data
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Return a dictionary with the bounding box data.

        EXAMPLES::

            sage: p = complex_plot(lambda z: z, (-1, 2), (-3, 4))
            sage: sorted(p.get_minmax_data().items())
            [('xmax', 2.0), ('xmin', -1.0), ('ymax', 4.0), ('ymin', -3.0)]
            sage: p = complex_plot(lambda z: z, (1, 2), (3, 4))
            sage: sorted(p.get_minmax_data().items())
            [('xmax', 2.0), ('xmin', 1.0), ('ymax', 4.0), ('ymin', 3.0)]
        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.x_range, self.y_range, dict=True)

    def _allowed_options(self):
        """
        TESTS::

            sage: isinstance(complex_plot(lambda z: z, (-1,1), (-1,1))[0]._allowed_options(), dict)
            True
        """
        return {'plot_points':'How many points to use for plotting precision',
                'interpolation':'What interpolation method to use'}

    def _repr_(self):
        """
        TESTS::

            sage: isinstance(complex_plot(lambda z: z, (-1,1), (-1,1))[0]._repr_(), str)
            True
        """
        return "ComplexPlot defined by a %s x %s data grid"%(self.x_count, self.y_count)

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: complex_plot(lambda x: x^2, (-5, 5), (-5, 5))
            Graphics object consisting of 1 graphics primitive
        """
        options = self.options()
        x0, x1 = float(self.x_range[0]), float(self.x_range[1])
        y0, y1 = float(self.y_range[0]), float(self.y_range[1])
        subplot.imshow(self.rgb_data, origin='lower', extent=(x0, x1, y0, y1),
                       interpolation=options['interpolation'])


@options(plot_points=100, interpolation='catrom')
def complex_plot(f, x_range, y_range, contoured=False, tiled=False, cmap=None,
                 contour_type='logarithmic', contour_base=None, dark_rate=0.5,
                 nphases=10, **options):
    r"""
    ``complex_plot`` takes a complex function of one variable,
    `f(z)` and plots output of the function over the specified
    ``x_range`` and ``y_range`` as demonstrated below. The magnitude of
    the output is indicated by the brightness and the argument is
    represented by the hue.

    By default, zero magnitude corresponds to black output, infinite
    magnitude corresponds to white output. The options ``contoured``,
    ``tiled``, and ``cmap`` affect the output.

    ``complex_plot(f, (xmin, xmax), (ymin, ymax), contoured, tiled, cmap, ...)``

    INPUT:

    - ``f`` -- a function of a single complex value `x + iy`

    - ``(xmin, xmax)`` -- 2-tuple, the range of ``x`` values

    - ``(ymin, ymax)`` -- 2-tuple, the range of ``y`` values

    - ``cmap`` -- ``None``, or the string name of a matplotlib colormap, or an
      instance of a matplotlib Colormap, or the special string ``'matplotlib'``
      (default: ``None``); if ``None``, then hues are chosen from a standard
      color wheel, cycling from red to yellow to blue. If ``matplotlib``, then
      hues are chosen from a preset matplotlib colormap.

    The following named parameter inputs can be used to add contours and adjust
    their distribution:

    - ``contoured`` -- boolean (default: ``False``); causes the magnitude
      to be indicated by logarithmically spaced 'contours'. The
      magnitude along one contour is either twice or half the magnitude
      along adjacent contours.

    - ``dark_rate`` -- a positive number (default: 0.5); affects how quickly
      magnitudes affect how light/dark the image is. When there are contours,
      this affects how visible each contour is. Large values (near 1.0) have
      very strong, immediate effects, while small values (near 0.0) have
      gradual effects.

    - ``tiled`` -- boolean (default: ``False``); causes the magnitude to
      be indicated by logarithmically spaced 'contours' as in
      ``contoured``, and in addition for there to be `10` evenly
      spaced phase contours.

    - ``nphases`` -- positive integer (default: 10); when ``tiled=True``,
      this is the number of divisions the phase is divided into

    - ``contour_type`` -- either ``'logarithmic'``, or ``'linear'`` (default:
      ``'logarithmic'``); causes added contours to be of given type when
      ``contoured=True``.

    - ``contour_base`` -- positive integer; when ``contour_type`` is
      ``'logarithmic'``, this sets logarithmic contours at multiples of
      ``contour_base`` apart. When ``contour_type`` is ``'linear'``, this sets
      contours at distances of ``contour_base`` apart. If ``None``, then a
      default is chosen depending on ``contour_type``.

    The following inputs may also be passed in as named parameters:

    - ``plot_points`` -- integer (default: 100); number of points to
      plot in each direction of the grid

    - ``interpolation`` -- string (default: ``'catrom'``); the interpolation
      method to use: ``'bilinear'``, ``'bicubic'``, ``'spline16'``,
      ``'spline36'``, ``'quadric'``, ``'gaussian'``, ``'sinc'``,
      ``'bessel'``, ``'mitchell'``, ``'lanczos'``, ``'catrom'``,
      ``'hermite'``, ``'hanning'``, ``'hamming'``, ``'kaiser'``

    Any additional parameters will be passed to ``show()``, as long as they're
    valid.

    .. NOTE::

        Matplotlib colormaps can be chosen or customized to cater to different
        types of vision. The colormaps 'cividis' and 'viridis' in matplotlib
        are designed to be perceptually uniform to a broader audience. The
        colormap 'turbo' is similar to the default but with more even contrast.
        See [NAR2018]_ for more information about colormap choice for
        scientific visualization.

    EXAMPLES:

    Here we plot a couple of simple functions::

        sage: complex_plot(sqrt(x), (-5, 5), (-5, 5))                                   # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(sqrt(x), (-5, 5), (-5, 5)))

    ::

        sage: complex_plot(sin(x), (-5, 5), (-5, 5))                                    # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(sin(x), (-5, 5), (-5, 5)))

    ::

        sage: complex_plot(log(x), (-10, 10), (-10, 10))                                # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(log(x), (-10, 10), (-10, 10)))

    ::

        sage: complex_plot(exp(x), (-10, 10), (-10, 10))                                # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(exp(x), (-10, 10), (-10, 10)))

    A plot with a different choice of colormap::

        sage: complex_plot(exp(x), (-10, 10), (-10, 10), cmap='viridis')                # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(exp(x), (-10, 10), (-10, 10), cmap='viridis'))

    A function with some nice zeros and a pole::

        sage: f(z) = z^5 + z - 1 + 1/z                                                  # needs sage.symbolic
        sage: complex_plot(f, (-3, 3), (-3, 3))                                         # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (-3, 3), (-3, 3)))

    The same function as above, but with contours. Contours render poorly with
    few plot points, so we use 300 here::

        sage: f(z) = z^5 + z - 1 + 1/z                                                  # needs sage.symbolic
        sage: complex_plot(f, (-3, 3), (-3, 3), plot_points=300, contoured=True)        # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (-3, 3), (-3, 3), plot_points=300, contoured=True))

    The same function as above, but tiled and with the *plasma* colormap::

        sage: f(z) = z^5 + z - 1 + 1/z                                                  # needs sage.symbolic
        sage: complex_plot(f, (-3, 3), (-3, 3),                                         # needs sage.symbolic
        ....:              plot_points=300, tiled=True, cmap='plasma')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (-3, 3), (-3, 3), plot_points=300, tiled=True, cmap='plasma'))

    When using ``tiled=True``, the number of phase subdivisions can be
    controlled by adjusting ``nphases``. We make the same plot with fewer
    tilings::

        sage: f(z) = z^5 + z - 1 + 1/z                                                  # needs sage.symbolic
        sage: complex_plot(f, (-3, 3), (-3, 3), plot_points=300,                        # needs sage.symbolic
        ....:              tiled=True, nphases=5, cmap='plasma')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (-3, 3), (-3, 3), plot_points=300, tiled=True, nphases=5, cmap='plasma'))

    It is also possible to use *linear* contours. We plot the same function
    above on an inset, setting contours to appear `1` apart::

        sage: f(z) = z^5 + z - 1 + 1/z                                                  # needs sage.symbolic
        sage: complex_plot(f, (0, 1), (0, 1), plot_points=300,                          # needs sage.symbolic
        ....:              contoured=True, contour_type='linear', contour_base=1)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (0, 1), (0, 1), plot_points=300, contoured=True, contour_type='linear', contour_base=1))

    Note that tightly spaced contours can lead to Moiré patterns and aliasing
    problems. For example::

        sage: f(z) = z^5 + z - 1 + 1/z                                                  # needs sage.symbolic
        sage: complex_plot(f, (-3, 3), (-3, 3), plot_points=300,                        # needs sage.symbolic
        ....:              contoured=True, contour_type='linear', contour_base=1)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (-3, 3), (-3, 3), plot_points=300, contoured=True, contour_type='linear', contour_base=1))

    When choosing colormaps, cyclic colormaps such as *twilight* or *hsv* might
    be considered more appropriate for showing changes in phase without sharp
    color contrasts::

        sage: f(z) = z^5 + z - 1 + 1/z                                                  # needs sage.symbolic
        sage: complex_plot(f, (-3, 3), (-3, 3), plot_points=300, cmap='twilight')       # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (-3, 3), (-3, 3), plot_points=300, cmap='twilight'))

    Passing *matplotlib* as the colormap gives a special colormap that is
    similar to the default::

        sage: f(z) = z^5 + z - 1 + 1/z                                                  # needs sage.symbolic
        sage: complex_plot(f, (-3, 3), (-3, 3),                                         # needs sage.symbolic
        ....:              plot_points=300, contoured=True, cmap='matplotlib')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (-3, 3), (-3, 3), plot_points=300, contoured=True, cmap='matplotlib'))

    Here is the identity, useful for seeing what values map to what colors::

        sage: complex_plot(lambda z: z, (-3, 3), (-3, 3))                               # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(lambda z: z, (-3, 3), (-3, 3)))

    The Riemann Zeta function::

        sage: complex_plot(zeta, (-30,30), (-30,30))                                    # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(zeta, (-30,30), (-30,30)))

    For advanced usage, it is possible to tweak many parameters. Increasing
    ``dark_rate`` will make regions become darker/lighter faster when there are no
    contours::

        sage: complex_plot(zeta, (-30, 30), (-30, 30), dark_rate=1.0)                   # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(zeta, (-30,30), (-30,30), dark_rate=1.0))

    Decreasing ``dark_rate`` has the opposite effect. When there are contours,
    adjust ``dark_rate`` affects how visible contours are. Compare::

        sage: complex_plot(zeta, (-1, 9), (10, 20), plot_points=200,    # long time, needs sage.symbolic
        ....:              contoured=True, cmap='twilight', dark_rate=0.2)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(zeta, (-1, 9), (10, 20), plot_points=200, contoured=True, cmap='twilight', dark_rate=0.2))

    and::

        sage: complex_plot(zeta, (-1, 9), (10, 20), plot_points=200,    # long time, needs sage.symbolic
        ....:              contoured=True, cmap='twilight', dark_rate=0.75)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(zeta, (-1, 9), (10, 20), plot_points=200, contoured=True, cmap='twilight', dark_rate=0.75))

    In practice, different values of ``dark_rate`` will work well with
    different colormaps.

    Extra options will get passed on to show(), as long as they are valid::

        sage: complex_plot(lambda z: z, (-3, 3), (-3, 3), figsize=[1,1])                # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    ::

        sage: complex_plot(lambda z: z, (-3, 3), (-3, 3)).show(figsize=[1,1])  # These are equivalent                   # needs sage.symbolic

    REFERENCES:

    Plotting complex functions with colormaps follows the strategy
    from [LD2021]_ and incorporates contour techniques described
    in [WegSem2010]_.

    TESTS:

    Test to make sure that using fast_callable functions works::

        sage: # needs sage.symbolic
        sage: f(x) = x^2
        sage: g = fast_callable(f, domain=CC, vars='x')
        sage: h = fast_callable(f, domain=CDF, vars='x')
        sage: P = complex_plot(f, (-10, 10), (-10, 10))
        sage: Q = complex_plot(g, (-10, 10), (-10, 10))
        sage: R = complex_plot(h, (-10, 10), (-10, 10))
        sage: S = complex_plot(exp(x)-sin(x), (-10, 10), (-10, 10))
        sage: P; Q; R; S
        Graphics object consisting of 1 graphics primitive
        Graphics object consisting of 1 graphics primitive
        Graphics object consisting of 1 graphics primitive
        Graphics object consisting of 1 graphics primitive

    Test to make sure symbolic functions still work without declaring
    a variable.  (We don't do this in practice because it doesn't use
    fast_callable, so it is much slower.)

    ::

        sage: complex_plot(sqrt, (-5, 5), (-5, 5))                                      # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive
    """
    import matplotlib as mpl
    import numpy as np
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid
    from sage.ext.fast_callable import fast_callable
    from sage.rings.complex_double import CDF

    try:
        f = fast_callable(f, domain=CDF, expect_one_var=True)
    except (AttributeError, TypeError, ValueError):
        pass

    cdef double x, y
    _, ranges = setup_for_eval_on_grid([], [x_range, y_range],
                                       options['plot_points'])
    x_range = ranges[0]
    y_range = ranges[1]
    cdef list z_values = []
    cdef list row
    for y in srange(*y_range, include_endpoint=True):
        row = []
        for x in srange(*x_range, include_endpoint=True):
            sig_check()
            row.append(f(new_CDF_element(x, y)))
        z_values.append(row)

    if cmap is None:
        # produce colors using the established default method
        rgbs = complex_to_rgb(
            z_values, contoured=contoured, tiled=tiled,
            contour_type=contour_type, contour_base=contour_base,
            dark_rate=dark_rate, nphases=nphases
        )
    else:
        # choose colors from colormap
        if isinstance(cmap, str):
            if cmap == 'matplotlib':
                domain = np.linspace(0, 1, 256)
                shifted_domain = np.roll(domain, 128)
                default_cmap = mpl.colors.LinearSegmentedColormap.from_list(
                    "sage_default", mpl.colormaps['hsv'](shifted_domain)
                )
                cmap = default_cmap
            else:
                cmap = mpl.colormaps[cmap]
        rgbs = complex_to_cmap_rgb(
            z_values, cmap=cmap, contoured=contoured, tiled=tiled,
            contour_type=contour_type, contour_base=contour_base,
            dark_rate=dark_rate, nphases=nphases
        )

    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options, ignore=['xmin', 'xmax']))
    g.add_primitive(ComplexPlot(rgbs, x_range[:2], y_range[:2], options))
    return g


def rgb_to_hls(rgb):
    r"""
    Convert array of rgb values (each in the range `[0, 1]`)
    to a numpy array of hls values (each in the range `[0, 1]`)

    INPUT:

    - ``rgb`` -- an `N \times 3` array of floats with values
      in the range `[0, 1]`; the rgb values at each point. (Note that the input
      can actually be of any dimension, such as `N \times M \times 3`, as long
      as the last dimension has length `3`).

    OUTPUT:

    An `N \times 3` Numpy array of floats in the range `[0, 1]`, with
    the same dimensions as the input array.

    .. SEEALSO::

        :func:`sage.plot.complex_plot.hls_to_rgb`

    EXAMPLES:

    We convert a row of floats and verify that we can convert back using
    ``hls_to_rgb``::

        sage: from sage.plot.complex_plot import rgb_to_hls, hls_to_rgb
        sage: rgb = [[0.2, 0.4, 0.5], [0.1, 0.3, 1.0]]
        sage: hls = rgb_to_hls(rgb)
        sage: hls   # abs tol 1e-4
        array([[0.55555556, 0.35      , 0.42857143],
               [0.62962963, 0.55      , 1.        ]])
        sage: hls_to_rgb(hls)  # abs tol 1e-4
        array([[0.2, 0.4, 0.5],
               [0.1, 0.3, 1. ]])

    Multidimensional inputs can be given as well::

        sage: multidim_arr = [[[0, 0.2, 0.4], [1, 1, 1]], [[0, 0, 0], [0.5, 0.6, 0.9]]]
        sage: rgb_to_hls(multidim_arr)   # abs tol 1e-4
        array([[[0.58333333, 0.2       , 1.        ],
                [0.        , 1.        , 0.        ]],
               [[0.        , 0.        , 0.        ],
                [0.625     , 0.7       , 0.66666667]]])
    """
    # Notation and algorithm corresponds to colorsys.rgb_to_hls
    import numpy as np
    rgb = np.asarray(rgb)
    if rgb.shape[-1] != 3:
        raise ValueError("Last dimension of input array must be 3; "
                         "shape {} was found.".format(rgb.shape))
    in_shape = rgb.shape
    rgb = np.asarray(rgb, dtype=np.dtype(float))
    rgb_max = rgb.max(-1)
    rgb_min = rgb.min(-1)
    l = (rgb_max + rgb_min)/2.0  # lightness

    hls = np.zeros_like(rgb)
    delta = np.ptp(rgb, -1)
    s = np.zeros_like(delta)

    ipos = delta > 0
    idx = (l <= 0.5) & ipos
    s[idx] = delta[idx] / (rgb_max[idx] + rgb_min[idx])

    idx = (l > 0.5) & ipos
    s[idx] = delta[idx] / (2.0 - rgb_max[idx] - rgb_min[idx])  # saturation

    # red is max
    idx = (rgb[..., 0] == rgb_max) & ipos
    hls[idx, 0] = (rgb[idx, 1] - rgb[idx, 2]) / delta[idx]

    # green is max
    idx = (rgb[..., 1] == rgb_max) & ipos
    hls[idx, 0] = 2.0 + (rgb[idx, 2] - rgb[idx, 0]) / delta[idx]

    # blue is max
    idx = (rgb[..., 2] == rgb_max) & ipos
    hls[idx, 0] = 4.0 + (rgb[idx, 0] - rgb[idx, 1]) / delta[idx]

    hls[..., 0] = (hls[..., 0] / 6.0) % 1.0
    hls[..., 1] = l
    hls[..., 2] = s

    return hls.reshape(in_shape)


def hls_to_rgb(hls):
    r"""
    Convert array of hls values (each in the range `[0, 1]`)
    to a numpy array of rgb values (each in the range `[0, 1]`)

    INPUT:

    - ``hls`` -- an `N \times 3` array of floats in the range `[0, 1]`; the hls
      values at each point. (Note that the input can actually be of any
      dimension, such as `N \times M \times 3`, as long as the last dimension
      has length `3`).

    OUTPUT:

    An `N \times 3`  Numpy array of floats in the range `[0, 1]`, with
    the same dimensions as the input array.

    .. SEEALSO::

        :func:`sage.plot.complex_plot.rgb_to_hls`

    EXAMPLES:

    We convert a row of floats and verify that we can convert back using
    ``rgb_to_hls``::

        sage: from sage.plot.complex_plot import rgb_to_hls, hls_to_rgb
        sage: hls = [[0.2, 0.4, 0.5], [0.1, 0.3, 1.0]]
        sage: rgb = hls_to_rgb(hls)
        sage: rgb  # abs tol 1e-4
        array([[0.52, 0.6 , 0.2 ],
               [0.6 , 0.36, 0.  ]])
        sage: rgb_to_hls(rgb)  # abs tol 1e-4
        array([[0.2, 0.4, 0.5],
               [0.1, 0.3, 1. ]])

    Multidimensional inputs can be given as well::

        sage: multidim_arr = [[[0, 0.2, 0.4], [0, 1, 0]], [[0, 0, 0], [0.5, 0.6, 0.9]]]
        sage: hls_to_rgb(multidim_arr)  # abs tol 1e-4
        array([[[0.28, 0.12, 0.12],
                [1.  , 1.  , 1.  ]],
               [[0.  , 0.  , 0.  ],
                [0.24, 0.96, 0.96]]])
    """
    import numpy as np
    hls = np.asarray(hls)
    if hls.shape[-1] != 3:
        raise ValueError("Last dimension of input array must be 3; "
                         "shape {} was found.".format(hls.shape))
    in_shape = hls.shape
    hls = np.array(
        hls, copy=False, dtype=np.dtype(float), ndmin=2
    )
    rgb = np.zeros_like(hls)

    szero_idx = hls[..., 2] == 0
    rgb[szero_idx, 0] = hls[szero_idx, 1]
    rgb[szero_idx, 1] = hls[szero_idx, 1]
    rgb[szero_idx, 2] = hls[szero_idx, 1]

    snonzero_idx = hls[..., 2] > 0
    m1 = np.zeros_like(snonzero_idx, dtype=float)
    m2 = np.zeros_like(snonzero_idx, dtype=float)

    idx = (hls[..., 1] <= 0.5)
    m2[idx] = hls[idx, 1] * (1.0 + hls[idx, 2])

    idx = (hls[..., 1] > 0.5)
    m2[idx] = hls[idx, 1] + hls[idx, 2] - hls[idx, 1] * hls[idx, 2]

    m1 = 2 * hls[..., 1] - m2

    rgb[snonzero_idx, 0] = _v(m1, m2, hls[..., 0] + 1.0 / 3.0)[snonzero_idx]
    rgb[snonzero_idx, 1] = _v(m1, m2, hls[..., 0])[snonzero_idx]
    rgb[snonzero_idx, 2] = _v(m1, m2, hls[..., 0] - 1.0 / 3.0)[snonzero_idx]

    return rgb.reshape(in_shape)


def _v(m1, m2, hue):
    """
    A helper function to convert from hls to rgb.

    This is a numpy version of colorsys._v.

    INPUT:

    - ``m1`` -- an array of floats with values in the range `[0, 1]`

    - ``m2`` -- an array of floats with values in the range `[0, 1]` with
      the same dimensions as ``m1``

    - ``hue`` -- an array of floats with values in the range `[0, 1]` with
      the same dimensions as ``m1``

    OUTPUT:

    A Numpy array of floats with values in the range `[0, 1]`. The dimensions
    are the same as the input arrays.

    EXAMPLES::

        sage: from sage.plot.complex_plot import _v
        sage: _v([0.1, 0.2], [0.1, 0.3], [0.25, 0.75])  # abs tol 1e-4
        array([0.1, 0.2])
    """
    import numpy as np
    m1 = np.asarray(m1, dtype=float)
    m2 = np.asarray(m2, dtype=float)
    hue = np.asarray(hue, dtype=float)

    if not m1.shape == m2.shape and m1.shape == hue.shape:
        raise ValueError(
            "Incompatible shapes given to _v. "
            "shapes {}, {}, {} given.".format(m1.shape, m2.shape, hue.shape)
        )

    out = np.zeros_like(m1, dtype=float)
    hue = hue % 1.0

    conditions = [hue < 1.0 / 6.0, hue < 1.0 / 2.0, hue < 2.0 / 3.0]
    out = np.select(conditions, [
        m1 + (m2 - m1) * hue * 6.0,
        m2,
        m1 + (m2 - m1) * (2.0 / 3.0 - hue) * 6.0
    ], default=m1)
    return out
