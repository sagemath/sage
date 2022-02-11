r"""
Spike Functions

AUTHORS:

- William Stein (2007-07): initial version

- Karl-Dieter Crisman (2009-09): adding documentation and doctests
"""
# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2009 Karl-Dieter Crisman <kcrisman@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import math

from sage.plot.all import line
from sage.modules.free_module_element import vector
from sage.rings.real_double import RDF


class SpikeFunction:
    """
    Base class for spike functions.

    INPUT:

    -  ``v`` - list of pairs (x, height)

    -  ``eps`` - parameter that determines approximation to a true spike

    OUTPUT:

    a function with spikes at each point ``x`` in ``v`` with the given height.

    EXAMPLES::

        sage: spike_function([(-3,4),(-1,1),(2,3)],0.001)
        A spike function with spikes at [-3.0, -1.0, 2.0]

    Putting the spikes too close together may delete some::

        sage: spike_function([(1,1),(1.01,4)],0.1)
        Some overlapping spikes have been deleted.
        You might want to use a smaller value for eps.
        A spike function with spikes at [1.0]

    Note this should normally be used indirectly via
    ``spike_function``, but one can use it directly::

        sage: from sage.functions.spike_function import SpikeFunction
        sage: S = SpikeFunction([(0,1),(1,2),(pi,-5)])
        sage: S
        A spike function with spikes at [0.0, 1.0, 3.141592653589793]
        sage: S.support
        [0.0, 1.0, 3.141592653589793]
    """
    def __init__(self, v, eps=0.0000001):
        """
        Initialize base class SpikeFunction.

        EXAMPLES::

            sage: S = spike_function([(-3,4),(-1,1),(2,3)],0.001); S
            A spike function with spikes at [-3.0, -1.0, 2.0]
            sage: S.height
            [4.0, 1.0, 3.0]
            sage: S.eps
            0.00100000000000000
        """
        if not v:
            v = [(0, 0)]
        v = sorted([(float(x[0]), float(x[1])) for x in v])
        notify = False

        for i in reversed(range(len(v)-1)):
            if v[i+1][0] - v[i][0] <= eps:
                notify = True
                del v[i+1]

        if notify:
            print("Some overlapping spikes have been deleted.")
            print("You might want to use a smaller value for eps.")

        self.v = v
        self.eps = eps
        self.support = [float(x[0]) for x in self.v]
        self.height = [float(x[1]) for x in self.v]

    def __repr__(self):
        """
        String representation of a spike function.

        EXAMPLES::

            sage: spike_function([(-3,4),(-1,1),(2,3)],0.001)
            A spike function with spikes at [-3.0, -1.0, 2.0]
        """
        return "A spike function with spikes at %s" % self.support

    def _eval(self, x):
        """
        Evaluate spike function.

        Note that when one calls the function within the tolerance,
        the return value is the full height at that point.

        EXAMPLES::

            sage: S = spike_function([(0,5)],eps=.001)
            sage: S(0)
            5.0
            sage: S(.1)
            0.0
            sage: S(.01)
            0.0
            sage: S(.001)
            5.0
        """
        eps = self.eps
        x = float(x)
        for i in range(len(self.support)):
            z = self.support[i]
            if z - eps <= x and x <= z + eps:
                return self.height[i], i
        return float(0), -1

    def __call__(self, x):
        """
        Called when spike function is used as callable function.

        EXAMPLES::

            sage: S = spike_function([(0,5)],eps=.001)
            sage: S(0)
            5.0
            sage: S(.1)
            0.0
            sage: S(.01)
            0.0
            sage: S(.001)
            5.0
        """
        return self._eval(x)[0]

    def plot_fft_abs(self, samples=2**12, xmin=None, xmax=None,  **kwds):
        """
        Plot of (absolute values of) Fast Fourier Transform of
        the spike function with given number of samples.

        EXAMPLES::

            sage: S = spike_function([(-3,4),(-1,1),(2,3)]); S
            A spike function with spikes at [-3.0, -1.0, 2.0]
            sage: P = S.plot_fft_abs(8)
            sage: p = P[0]; p.ydata  # abs tol 1e-8
            [5.0, 5.0, 3.367958691924177, 3.367958691924177, 4.123105625617661,
             4.123105625617661, 4.759921664218055, 4.759921664218055]
        """
        w = self.vector(samples=samples, xmin=xmin, xmax=xmax)
        z = w.fft()
        k = vector(RDF, [abs(z[i]) for i in range(len(z)//2)])
        return k.plot(xmin=0, xmax=1, **kwds)

    def plot_fft_arg(self, samples=2**12, xmin=None, xmax=None,  **kwds):
        """
        Plot of (absolute values of) Fast Fourier Transform of
        the spike function with given number of samples.

        EXAMPLES::

            sage: S = spike_function([(-3,4),(-1,1),(2,3)]); S
            A spike function with spikes at [-3.0, -1.0, 2.0]
            sage: P = S.plot_fft_arg(8)
            sage: p = P[0]; p.ydata  # abs tol 1e-8
            [0.0, 0.0, -0.211524990023434, -0.211524990023434,
             0.244978663126864, 0.244978663126864, -0.149106180027477,
             -0.149106180027477]
        """
        w = self.vector(samples=samples, xmin=xmin, xmax=xmax)
        z = w.fft()
        k = vector(RDF, [(z[i]).arg() for i in range(len(z)//2)])
        return k.plot(xmin=0, xmax=1, **kwds)

    def vector(self, samples=2**16, xmin=None, xmax=None):
        """
        Create a sampling vector of the spike function in question.

        EXAMPLES::

            sage: S = spike_function([(-3,4),(-1,1),(2,3)],0.001); S
            A spike function with spikes at [-3.0, -1.0, 2.0]
            sage: S.vector(16)
            (4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0)
        """
        v = vector(RDF, samples)  # creates vector of zeros of length 2^16
        xmin, xmax = self._ranges(xmin, xmax)
        delta = (xmax - xmin) / samples
        w = int(math.ceil(self.eps / delta))
        for i in range(len(self.support)):
            x = self.support[i]
            if x > xmax:
                break
            h = self.height[i]
            j = int((x - xmin) / delta)
            for k in range(j, min(samples, j + w)):
                v[k] = h
        return v

    def _ranges(self, xmin, xmax):
        """
        Quickly find appropriate plotting interval.

        EXAMPLES::

            sage: S = spike_function([(-1,1),(1,40)])
            sage: S._ranges(None,None)
            (-1.0, 1.0)
        """
        width = (self.support[-1] + self.support[0])/float(2)
        if xmin is None:
            xmin = self.support[0] - width/float(5)
        if xmax is None:
            xmax = self.support[-1] + width/float(5)
        if xmax <= xmin:
            xmax = xmin + 1
        return xmin, xmax

    def plot(self, xmin=None, xmax=None, **kwds):
        """
        Special fast plot method for spike functions.

        EXAMPLES::

            sage: S = spike_function([(-1,1),(1,40)])
            sage: P = plot(S)
            sage: P[0]
            Line defined by 8 points
        """
        v = []
        xmin, xmax = self._ranges(xmin, xmax)
        x = xmin
        eps = self.eps
        while x < xmax:
            y, i = self._eval(x)
            v.append((x, y))
            if i != -1:
                x0 = self.support[i] + eps
                v.extend([(x0, y), (x0, 0)])
                if i+1 < len(self.support):
                    x = self.support[i + 1] - eps
                    v.append((x, 0))
                else:
                    x = xmax
                    v.append((xmax, 0))
            else:
                new_x = None
                for j in range(len(self.support)):
                    if self.support[j] - eps > x:
                        new_x = self.support[j] - eps
                        break
                if new_x is None:
                    new_x = xmax
                v.append((new_x, 0))
                x = new_x
        L = line(v, **kwds)
        L.xmin(xmin-1)
        L.xmax(xmax)
        return L


spike_function = SpikeFunction
