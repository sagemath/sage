r"""
Discrete Fourier Transforms

This file contains functions useful for computing discrete Fourier
transforms and probability distribution functions for discrete random
variables for sequences of elements of `\QQ` or `\CC`, indexed by a
``range(N)``, `\ZZ / N \ZZ`, an abelian group, the conjugacy classes
of a permutation group, or the conjugacy classes of a matrix group.

This file implements:

- :meth:`__eq__`

- :meth:`__mul__` (for right multiplication by a scalar)

- plotting, printing -- :meth:`IndexedSequence.plot`,
  :meth:`IndexedSequence.plot_histogram`, :meth:`_repr_`, :meth:`__str__`

- dft --  computes the discrete Fourier transform for the following cases:

  * a sequence (over `\QQ` or :class:`CyclotomicField`) indexed by ``range(N)``
    or `\ZZ / N \ZZ`
  * a sequence (as above) indexed by a finite abelian group
  * a sequence (as above) indexed by a complete set of representatives of
    the conjugacy classes of a finite permutation group
  * a sequence (as above) indexed by a complete set of representatives of
    the conjugacy classes of a finite matrix group

- idft --  computes the discrete Fourier transform for the following cases:

  * a sequence (over `\QQ` or CyclotomicField) indexed by ``range(N)`` or
    `\ZZ / N \ZZ`

- dct, dst  (for discrete Fourier/Cosine/Sine transform)

- convolution (in :meth:`IndexedSequence.convolution` and
  :meth:`IndexedSequence.convolution_periodic`)

- fft, ifft -- (fast Fourier transforms) wrapping GSL's
  ``gsl_fft_complex_forward()``, ``gsl_fft_complex_inverse()``,
  using William Stein's :func:`FastFourierTransform`

- dwt, idwt -- (fast wavelet transforms) wrapping GSL's ``gsl_dwt_forward()``,
  ``gsl_dwt_backward()`` using Joshua Kantor's :func:`WaveletTransform` class.
  Allows for wavelets of type:

  * "haar"
  * "daubechies"
  * "daubechies_centered"
  * "haar_centered"
  * "bspline"
  * "bspline_centered"


.. TODO::

    - "filtered" DFTs
    - more idfts
    - more examples for probability, stats, theory of FTs

AUTHORS:

- David Joyner (2006-10)

- William Stein (2006-11) -- fix many bugs
"""

##########################################################################
#  Copyright (C) 2006 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL):
#
#                  http://www.gnu.org/licenses/
##########################################################################

from sage.misc.superseded import deprecation
deprecation(9084, "the module sage.gsl.dft has moved to sage.calculus.transforms.dft")

from sage.calculus.transforms.dft import *
