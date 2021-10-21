"""
Basic Statistics

This file contains basic descriptive functions. Included are the mean,
median, mode, moving average, standard deviation, and the variance.
When calling a function on data, there are checks for functions already
defined for that data type.

The ``mean`` function returns the arithmetic mean (the sum of all the members
of a list, divided by the number of members). Further revisions may include
the geometric and harmonic mean. The ``median`` function returns the number
separating the higher half of a sample from the lower half. The ``mode``
returns the most common occurring member of a sample, plus the number of times
it occurs. If entries occur equally common, the smallest of a list of the most
common  entries is returned. The ``moving_average`` is a finite impulse
response filter, creating a series of averages using a user-defined number of
subsets of the full data set. The ``std`` and the ``variance`` return a
measurement of how far data points tend to be from the arithmetic mean.

Functions are available in the namespace ``stats``, i.e. you can use them by
typing ``stats.mean``, ``stats.median``, etc.

REMARK: If all the data you are working with are floating point
numbers, you may find ``stats.TimeSeries`` helpful, since it is
extremely fast and offers many of the same descriptive statistics as
in the module.

AUTHOR:

- Andrew Hou (11/06/2009)

"""
# ***********************************************************************
#          Copyright (C) 2009, Andrew Hou <amhou@uw.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#            The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
######################################################################

from sage.rings.integer_ring import ZZ
from sage.symbolic.constants import NaN
from sage.misc.functional import sqrt
from sage.misc.superseded import deprecation


def mean(v):
    """
    Return the mean of the elements of `v`.

    We define the mean of the empty list to be the (symbolic) NaN,
    following the convention of MATLAB, Scipy, and R.

    This function is deprecated.  Use ``numpy.mean`` or ``numpy.nanmean``
    instead.

    INPUT:

    - `v` -- a list of numbers

    OUTPUT:

    - a number

    EXAMPLES::

        sage: mean([pi, e])
        doctest:warning...
        DeprecationWarning: sage.stats.basic_stats.mean is deprecated; use numpy.mean or numpy.nanmean instead
        See https://trac.sagemath.org/29662 for details.
        1/2*pi + 1/2*e
        sage: mean([])
        NaN
        sage: mean([I, sqrt(2), 3/5])
        1/3*sqrt(2) + 1/3*I + 1/5
        sage: mean([RIF(1.0103,1.0103), RIF(2)])
        1.5051500000000000?
        sage: mean(range(4))
        3/2
        sage: v = stats.TimeSeries([1..100])
        sage: mean(v)
        50.5
    """
    deprecation(29662, 'sage.stats.basic_stats.mean is deprecated; use numpy.mean or numpy.nanmean instead')
    if hasattr(v, 'mean'):
        return v.mean()
    if not v:
        return NaN
    s = sum(v)
    if isinstance(s, int):
        # python integers are stupid.
        return s / ZZ(len(v))
    return s / len(v)


def mode(v):
    """
    Return the mode of `v`.

    The mode is the list of the most frequently occurring
    elements in `v`. If `n` is the most times that any element occurs
    in `v`, then the mode is the list of elements of `v` that
    occur `n` times. The list is sorted if possible.

    This function is deprecated.  Use ``scipy.stats.mode`` or
    ``statistics.mode`` instead.

    .. NOTE::

        The elements of `v` must be hashable.

    INPUT:

    - `v` -- a list

    OUTPUT:

    - a list (sorted if possible)

    EXAMPLES::

        sage: v = [1,2,4,1,6,2,6,7,1]
        sage: mode(v)
        doctest:warning...
        DeprecationWarning: sage.stats.basic_stats.mode is deprecated; use scipy.stats.mode or statistics.mode instead
        See https://trac.sagemath.org/29662 for details.
        [1]
        sage: v.count(1)
        3
        sage: mode([])
        []

        sage: mode([1,2,3,4,5])
        [1, 2, 3, 4, 5]
        sage: mode([3,1,2,1,2,3])
        [1, 2, 3]
        sage: mode([0, 2, 7, 7, 13, 20, 2, 13])
        [2, 7, 13]

        sage: mode(['sage', 'four', 'I', 'three', 'sage', 'pi'])
        ['sage']

        sage: class MyClass:
        ....:   def mode(self):
        ....:       return [1]
        sage: stats.mode(MyClass())
        [1]
    """
    deprecation(29662, 'sage.stats.basic_stats.mode is deprecated; use scipy.stats.mode or statistics.mode instead')

    if hasattr(v, 'mode'):
        return v.mode()

    if not v:
        return v

    freq = {}
    for i in v:
        if i in freq:
            freq[i] += 1
        else:
            freq[i] = 1

    n = max(freq.values())
    try:
        return sorted(u for u, f in freq.items() if f == n)
    except TypeError:
        return [u for u, f in freq.items() if f == n]


def std(v, bias=False):
    """
    Return the standard deviation of the elements of `v`.

    We define the standard deviation of the empty list to be NaN,
    following the convention of MATLAB, Scipy, and R.

    This function is deprecated.  Use ``numpy.std`` or ``numpy.nanstd``
    instead.

    INPUT:

    - `v` -- a list of numbers

    - ``bias`` -- bool (default: False); if False, divide by
                  len(v) - 1 instead of len(v)
                  to give a less biased estimator (sample) for the
                  standard deviation.

    OUTPUT:

    - a number

    EXAMPLES::

        sage: std([1..6], bias=True)
        doctest:warning...
        DeprecationWarning: sage.stats.basic_stats.std is deprecated; use numpy.std or numpy.nanstd instead
        See https://trac.sagemath.org/29662 for details.
        doctest:warning...
        DeprecationWarning: sage.stats.basic_stats.variance is deprecated; use numpy.var or numpy.nanvar instead
        See https://trac.sagemath.org/29662 for details.
        doctest:warning...
        DeprecationWarning: sage.stats.basic_stats.mean is deprecated; use numpy.mean or numpy.nanmean instead
        See https://trac.sagemath.org/29662 for details.
        1/2*sqrt(35/3)
        sage: std([1..6], bias=False)
        sqrt(7/2)
        sage: std([e, pi])
        sqrt(1/2)*abs(pi - e)
        sage: std([])
        NaN
        sage: std([I, sqrt(2), 3/5])
        1/15*sqrt(1/2)*sqrt((10*sqrt(2) - 5*I - 3)^2
        + (5*sqrt(2) - 10*I + 3)^2 + (5*sqrt(2) + 5*I - 6)^2)
        sage: std([RIF(1.0103, 1.0103), RIF(2)])
        0.6998235813403261?
        sage: import numpy
        sage: x = numpy.array([1,2,3,4,5])
        sage: std(x, bias=False)
        1.5811388300841898
        sage: x = stats.TimeSeries([1..100])
        sage: std(x)
        29.011491975882016

    TESTS::

        sage: data = [random() for i in [1 .. 20]]
        sage: std(data)  # random
        0.29487771726609185
    """
    deprecation(29662, 'sage.stats.basic_stats.std is deprecated; use numpy.std or numpy.nanstd instead')

    # NOTE: in R bias = False by default, and in Scipy bias=True by
    # default, and R is more popular.

    if hasattr(v, 'standard_deviation'):
        return v.standard_deviation(bias=bias)

    import numpy

    if isinstance(v, numpy.ndarray):
        # accounts for numpy arrays
        if bias:
            return v.std()
        else:
            return v.std(ddof=1)

    if not v:
        # standard deviation of empty set defined as NaN
        return NaN

    return sqrt(variance(v, bias=bias))


def variance(v, bias=False):
    """
    Return the variance of the elements of `v`.

    We define the variance of the empty list to be NaN,
    following the convention of MATLAB, Scipy, and R.

    This function is deprecated.  Use ``numpy.var`` or ``numpy.nanvar``
    instead.

    INPUT:

    - `v` -- a list of numbers

    - ``bias`` -- bool (default: False); if False, divide by
                  len(v) - 1 instead of len(v)
                  to give a less biased estimator (sample) for the
                  standard deviation.

    OUTPUT:

    - a number

    EXAMPLES::

        sage: variance([1..6])
        doctest:warning...
        DeprecationWarning: sage.stats.basic_stats.variance is deprecated; use numpy.var or numpy.nanvar instead
        See https://trac.sagemath.org/29662 for details.
        7/2
        sage: variance([1..6], bias=True)
        35/12
        sage: variance([e, pi])
        1/2*(pi - e)^2
        sage: variance([])
        NaN
        sage: variance([I, sqrt(2), 3/5])
        1/450*(10*sqrt(2) - 5*I - 3)^2 + 1/450*(5*sqrt(2) - 10*I + 3)^2
        + 1/450*(5*sqrt(2) + 5*I - 6)^2
        sage: variance([RIF(1.0103, 1.0103), RIF(2)])
        0.4897530450000000?
        sage: import numpy
        sage: x = numpy.array([1,2,3,4,5])
        sage: variance(x, bias=False)
        2.5
        sage: x = stats.TimeSeries([1..100])
        sage: variance(x)
        841.6666666666666
        sage: variance(x, bias=True)
        833.25
        sage: class MyClass:
        ....:   def variance(self, bias = False):
        ....:      return 1
        sage: stats.variance(MyClass())
        1
        sage: class SillyPythonList:
        ....:   def __init__(self):
        ....:       self.__list = [2, 4]
        ....:   def __len__(self):
        ....:       return len(self.__list)
        ....:   def __iter__(self):
        ....:       return self.__list.__iter__()
        ....:   def mean(self):
        ....:       return 3
        sage: R = SillyPythonList()
        sage: variance(R)
        2
        sage: variance(R, bias=True)
        1

    TESTS:

    The performance issue from :trac:`10019` is solved::

        sage: variance([1] * 2^18)
        0
    """
    deprecation(29662, 'sage.stats.basic_stats.variance is deprecated; use numpy.var or numpy.nanvar instead')

    if hasattr(v, 'variance'):
        return v.variance(bias=bias)
    import numpy

    x = 0
    if isinstance(v, numpy.ndarray):
        # accounts for numpy arrays
        if bias:
            return v.var()
        else:
            return v.var(ddof=1)
    if not v:
        # variance of empty set defined as NaN
        return NaN

    mu = mean(v)
    for vi in v:
        x += (vi - mu)**2
    if bias:
        # population variance
        if isinstance(x, int):
            return x / ZZ(len(v))
        return x / len(v)
    else:
        # sample variance
        if isinstance(x, int):
            return x / ZZ(len(v)-1)
        return x / (len(v)-1)


def median(v):
    """
    Return the median (middle value) of the elements of `v`

    If `v` is empty, we define the median to be NaN, which is
    consistent with NumPy (note that R returns NULL).
    If `v` is comprised of strings, TypeError occurs.
    For elements other than numbers, the median is a result of ``sorted()``.

    This function is deprecated.  Use ``numpy.median`` or ``numpy.nanmedian``
    instead.

    INPUT:

    - `v` -- a list

    OUTPUT:

    - median element of `v`

    EXAMPLES::

        sage: median([1,2,3,4,5])
        doctest:warning...
        DeprecationWarning: sage.stats.basic_stats.median is deprecated; use numpy.median or numpy.nanmedian instead
        See https://trac.sagemath.org/29662 for details.
        3
        sage: median([e, pi])
        1/2*pi + 1/2*e
        sage: median(['sage', 'linux', 'python'])
        'python'
        sage: median([])
        NaN
        sage: class MyClass:
        ....:    def median(self):
        ....:       return 1
        sage: stats.median(MyClass())
        1
    """
    deprecation(29662, 'sage.stats.basic_stats.median is deprecated; use numpy.median or numpy.nanmedian instead')

    if hasattr(v, 'median'):
        return v.median()

    if not v:
        # Median of empty set defined as NaN
        return NaN
    values = sorted(v)
    if len(values) % 2:
        return values[((len(values))+1)//2-1]
    else:
        lower = values[(len(values)+1)//2-1]
        upper = values[len(values)//2]
        return (lower + upper) / ZZ(2)


def moving_average(v, n):
    """
    Return the moving average of a list `v`.

    The moving average of a list is often used to smooth out noisy data.

    If `v` is empty, we define the entries of the moving average to be NaN.

    This method is deprecated.  Use ``pandas.Series.rolling`` instead.

    INPUT:

    - `v` -- a list

    - `n` -- the number of values used in computing each average.

    OUTPUT:

    - a list of length ``len(v)-n+1``, since we do not fabric any values

    EXAMPLES::

        sage: moving_average([1..10], 1)
        doctest:warning...
        DeprecationWarning: sage.stats.basic_stats.moving_average is deprecated; use pandas.Series.rolling instead
        See https://trac.sagemath.org/29662 for details.
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        sage: moving_average([1..10], 4)
        [5/2, 7/2, 9/2, 11/2, 13/2, 15/2, 17/2]
        sage: moving_average([], 1)
        []
        sage: moving_average([pi, e, I, sqrt(2), 3/5], 2)
        [1/2*pi + 1/2*e, 1/2*e + 1/2*I, 1/2*sqrt(2) + 1/2*I,
         1/2*sqrt(2) + 3/10]

    We check if the input is a time series, and if so use the
    optimized ``simple_moving_average`` method, but with (slightly
    different) meaning as defined above (the point is that the
    ``simple_moving_average`` on time series returns `n` values::

        sage: a = stats.TimeSeries([1..10])
        sage: stats.moving_average(a, 3)
        [2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000]
        sage: stats.moving_average(list(a), 3)
        [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

    """
    deprecation(29662, 'sage.stats.basic_stats.moving_average is deprecated; use pandas.Series.rolling instead')

    if not v:
        return v
    from .time_series import TimeSeries
    if isinstance(v, TimeSeries):
        return v.simple_moving_average(n)[n - 1:]
    n = int(n)
    if n <= 0:
        raise ValueError("n must be positive")
    nn = ZZ(n)
    s = sum(v[:n])
    ans = [s / nn]
    for i in range(n, len(v)):
        # add in the i-th value in v to our running sum,
        # and remove the value n places back.
        s += v[i] - v[i - n]
        ans.append(s / nn)
    return ans
