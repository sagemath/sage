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

- :meth:`dft` --  computes the discrete Fourier transform for the following cases:

  * a sequence (over `\QQ` or :class:`CyclotomicField`) indexed by ``range(N)``
    or `\ZZ / N \ZZ`
  * a sequence (as above) indexed by a finite abelian group
  * a sequence (as above) indexed by a complete set of representatives of
    the conjugacy classes of a finite permutation group
  * a sequence (as above) indexed by a complete set of representatives of
    the conjugacy classes of a finite matrix group

- :meth:`idft` --  computes the discrete Fourier transform for the following cases:

  * a sequence (over `\QQ` or CyclotomicField) indexed by ``range(N)`` or
    `\ZZ / N \ZZ`

- :meth:`dct`, :meth:`dst`  (for discrete Fourier/Cosine/Sine transform)

- convolution (in :meth:`IndexedSequence.convolution` and
  :meth:`IndexedSequence.convolution_periodic`)

- :meth:`fft`, :meth:`ifft` -- (fast Fourier transforms) wrapping GSL's
  ``gsl_fft_complex_forward()``, ``gsl_fft_complex_inverse()``,
  using William Stein's :func:`FastFourierTransform`

- :meth:`dwt`, :meth:`idwt` -- (fast wavelet transforms) wrapping GSL's ``gsl_dwt_forward()``,
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
#                  https://www.gnu.org/licenses/
##########################################################################
from sage.functions.trig import sin, cos
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage.structure.sequence import Sequence

lazy_import("sage.calculus.transforms.dwt", "WaveletTransform")
lazy_import("sage.calculus.transforms.fft", "FastFourierTransform")
lazy_import("sage.groups.abelian_gps.abelian_group", "AbelianGroup")
lazy_import("sage.groups.perm_gps.permgroup_element", "PermutationGroupElement")
lazy_import("sage.plot.all", ["polygon", "line", "text"])
lazy_import("sage.rings.number_field.number_field", "CyclotomicField")


class IndexedSequence(SageObject):
    """
    An indexed sequence.

    INPUT:

    - ``L`` -- list

    - ``index_object`` -- must be a Sage object with an ``__iter__`` method
      containing the same number of elements as ``self``, which is a
      list of elements taken from a field
    """
    def __init__(self, L, index_object):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: J = list(range(10))
            sage: A = [1/10 for j in J]
            sage: s = IndexedSequence(A,J)
            sage: s
            Indexed sequence: [1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10]
                indexed by [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: s.dict()
            {0: 1/10,
             1: 1/10,
             2: 1/10,
             3: 1/10,
             4: 1/10,
             5: 1/10,
             6: 1/10,
             7: 1/10,
             8: 1/10,
             9: 1/10}
            sage: s.list()
            [1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10]
            sage: s.index_object()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: s.base_ring()
            Rational Field
        """
        try:
            ind = index_object.list()
        except AttributeError:
            ind = list(index_object)
        self._index_object = index_object
        self._list = Sequence(L)
        self._base_ring = self._list.universe()
        dict = {}
        for i in range(len(ind)):
            dict[ind[i]] = L[i]
        self._dict = dict

    def dict(self):
        """
        Return a python dict of ``self`` where the keys are elements in the
        indexing set.

        EXAMPLES::

            sage: J = list(range(10))
            sage: A = [1/10 for j in J]
            sage: s = IndexedSequence(A,J)
            sage: s.dict()
            {0: 1/10, 1: 1/10, 2: 1/10, 3: 1/10, 4: 1/10, 5: 1/10, 6: 1/10, 7: 1/10, 8: 1/10, 9: 1/10}
        """
        return self._dict

    def list(self):
        """
        Return the list of ``self``.

        EXAMPLES::

            sage: J = list(range(10))
            sage: A = [1/10 for j in J]
            sage: s = IndexedSequence(A,J)
            sage: s.list()
            [1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10]
        """
        return self._list

    def base_ring(self):
        r"""
        This just returns the common parent `R` of the `N` list
        elements. In some applications (say, when computing the
        discrete Fourier transform, dft), it is more accurate to think
        of the base_ring as the group ring `\QQ(\zeta_N)[R]`.

        EXAMPLES::

            sage: J = list(range(10))
            sage: A = [1/10 for j in J]
            sage: s = IndexedSequence(A,J)
            sage: s.base_ring()
            Rational Field
        """
        return self._base_ring

    def index_object(self):
        """
        Return the indexing object.

        EXAMPLES::

            sage: J = list(range(10))
            sage: A = [1/10 for j in J]
            sage: s = IndexedSequence(A,J)
            sage: s.index_object()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        return self._index_object

    def _repr_(self):
        """
        Implement print method.

        EXAMPLES::

            sage: A = [ZZ(i) for i in range(3)]
            sage: I = list(range(3))
            sage: s = IndexedSequence(A,I)
            sage: s
            Indexed sequence: [0, 1, 2]
             indexed by [0, 1, 2]
            sage: print(s)
            Indexed sequence: [0, 1, 2]
             indexed by [0, 1, 2]
            sage: I = GF(3)
            sage: A = [i^2 for i in I]
            sage: s = IndexedSequence(A,I); s
            Indexed sequence: [0, 1, 1]
             indexed by Finite Field of size 3
        """
        return "Indexed sequence: " + str(self.list()) + "\n    indexed by " + str(self.index_object())

    def plot_histogram(self, clr=(0, 0, 1), eps=0.4):
        r"""
        Plot the histogram plot of the sequence.

        The sequence is assumed to be real or from a finite field,
        with a real indexing set ``I`` coercible into `\RR`.

        Options are ``clr``, which is an RGB value, and ``eps``, which
        is the spacing between the bars.

        EXAMPLES::

            sage: J = list(range(3))
            sage: A = [ZZ(i^2)+1 for i in J]
            sage: s = IndexedSequence(A,J)
            sage: P = s.plot_histogram()                                                # needs sage.plot
            sage: show(P)                       # not tested                            # needs sage.plot
        """
        from sage.rings.real_mpfr import RR
        # elements must be coercible into RR
        I = self.index_object()
        N = len(I)
        S = self.list()
        P = [polygon([[RR(I[i]) - eps, 0],
                      [RR(I[i]) - eps, RR(S[i])],
                      [RR(I[i]) + eps, RR(S[i])],
                      [RR(I[i]) + eps, 0],
                      [RR(I[i]), 0]], rgbcolor=clr) for i in range(N)]
        T = [text(str(I[i]), (RR(I[i]), -0.8), fontsize=15, rgbcolor=(1, 0, 0))
             for i in range(N)]
        return sum(P) + sum(T)

    def plot(self):
        """
        Plot the points of the sequence.

        Elements of the sequence are assumed to be real or from a
        finite field, with a real indexing set ``I = range(len(self))``.

        EXAMPLES::

            sage: I = list(range(3))
            sage: A = [ZZ(i^2)+1 for i in I]
            sage: s = IndexedSequence(A,I)
            sage: P = s.plot()                                                          # needs sage.plot
            sage: show(P)                       # not tested                            # needs sage.plot
        """
        from sage.rings.real_mpfr import RR
        # elements must be coercible into RR
        I = self.index_object()
        S = self.list()
        return line([[RR(I[i]), RR(S[i])] for i in range(len(I) - 1)])

    def dft(self, chi=None):
        r"""
        A discrete Fourier transform "over `\QQ`" using exact
        `N`-th roots of unity.

        EXAMPLES::

            sage: J = list(range(6))
            sage: A = [ZZ(1) for i in J]
            sage: s = IndexedSequence(A,J)
            sage: s.dft(lambda x: x^2)                                                  # needs sage.rings.number_field
            Indexed sequence: [6, 0, 0, 6, 0, 0]
             indexed by [0, 1, 2, 3, 4, 5]
            sage: s.dft()                                                               # needs sage.rings.number_field
            Indexed sequence: [6, 0, 0, 0, 0, 0]
             indexed by [0, 1, 2, 3, 4, 5]

            sage: # needs sage.combinat sage.groups
            sage: G = SymmetricGroup(3)
            sage: J = G.conjugacy_classes_representatives()
            sage: s = IndexedSequence([1,2,3], J)  # 1,2,3 are the values of a class fcn on G
            sage: s.dft()   # the "scalar-valued Fourier transform" of this class fcn
            Indexed sequence: [8, 2, 2]
             indexed by [(), (1,2), (1,2,3)]

            sage: # needs sage.rings.number_field
            sage: J = AbelianGroup(2, [2,3], names='ab')
            sage: s = IndexedSequence([1,2,3,4,5,6], J)
            sage: s.dft()   # the precision of output is somewhat random and architecture dependent.
            Indexed sequence: [21.0000000000000,
                               -2.99999999999997 - 1.73205080756885*I,
                               -2.99999999999999 + 1.73205080756888*I,
                               -9.00000000000000 + 0.0000000000000485744257349999*I,
                               -0.00000000000000976996261670137 - 0.0000000000000159872115546022*I,
                               -0.00000000000000621724893790087 - 0.0000000000000106581410364015*I]
             indexed by Multiplicative Abelian group isomorphic to C2 x C3

            sage: # needs sage.groups sage.rings.number_field
            sage: J = CyclicPermutationGroup(6)
            sage: s = IndexedSequence([1,2,3,4,5,6], J)
            sage: s.dft()   # the precision of output is somewhat random and architecture dependent.
            Indexed sequence: [21.0000000000000,
                               -2.99999999999997 - 1.73205080756885*I,
                               -2.99999999999999 + 1.73205080756888*I,
                               -9.00000000000000 + 0.0000000000000485744257349999*I,
                               -0.00000000000000976996261670137 - 0.0000000000000159872115546022*I,
                               -0.00000000000000621724893790087 - 0.0000000000000106581410364015*I]
             indexed by Cyclic group of order 6 as a permutation group

            sage: # needs sage.rings.number_field
            sage: p = 7; J = list(range(p)); A = [kronecker_symbol(j,p) for j in J]
            sage: s = IndexedSequence(A, J)
            sage: Fs = s.dft()
            sage: c = Fs.list()[1]; [x/c for x in Fs.list()]; s.list()
            [0, 1, 1, -1, 1, -1, -1]
            [0, 1, 1, -1, 1, -1, -1]

        The DFT of the values of the quadratic residue symbol is itself, up to
        a constant factor (denoted c on the last line above).

        .. TODO::

            Read the parent of the elements of S; if `\QQ` or `\CC` leave as
            is; if AbelianGroup, use abelian_group_dual; if some other
            implemented Group (permutation, matrix), call .characters()
            and test if the index list is the set of conjugacy classes.
        """
        if chi is None:
            chi = lambda x: x
        J = self.index_object()   # index set of length N
        N = len(J)
        S = self.list()
        F = self.base_ring()   # elements must be coercible into QQ(zeta_N)
        if J[0] not in ZZ:
            G = J[0].parent()  # if J is not a range it is a group G
        if J[0] in ZZ and F.base_ring().fraction_field() == QQ:
            # assumes J is range(N)
            zeta = CyclotomicField(N).gen()
            FT = [sum([S[i] * chi(zeta**(i * j)) for i in J]) for j in J]
        elif (J[0] not in ZZ) and G.is_abelian() and F == ZZ or (F.is_field() and F.base_ring() == QQ):
            if isinstance(J[0], PermutationGroupElement):
                # J is a CyclicPermGp
                n = G.order()
                a = list(n.factor())
                invs = [x[0]**x[1] for x in a]
                G = AbelianGroup(len(a), invs)
            # assumes J is AbelianGroup(...)
            Gd = G.dual_group()
            FT = [sum([S[i] * chid(G.list()[i]) for i in range(N)])
                  for chid in Gd]
        elif (J[0] not in ZZ) and G.is_finite() and F == ZZ or (F.is_field() and F.base_ring() == QQ):
            # assumes J is the list of conj class representatives of a
            # PermutationGroup(...) or Matrixgroup(...)
            chi = G.character_table()
            FT = [sum([S[i] * chi[i, j] for i in range(N)]) for j in range(N)]
        else:
            raise ValueError(f"list elements must be in QQ(zeta_{N})")
        return IndexedSequence(FT, J)

    def idft(self):
        r"""
        A discrete inverse Fourier transform. Only works over `\QQ`.

        EXAMPLES::

            sage: J = list(range(5))
            sage: A = [ZZ(1) for i in J]
            sage: s = IndexedSequence(A,J)
            sage: fs = s.dft(); fs                                                      # needs sage.rings.number_field
            Indexed sequence: [5, 0, 0, 0, 0]
                indexed by [0, 1, 2, 3, 4]
            sage: it = fs.idft(); it                                                    # needs sage.rings.number_field
            Indexed sequence: [1, 1, 1, 1, 1]
                indexed by [0, 1, 2, 3, 4]
            sage: it == s                                                               # needs sage.rings.number_field
            True
        """
        F = self.base_ring()   # elements must be coercible into QQ(zeta_N)
        J = self.index_object()   # must be = range(N)
        N = len(J)
        S = self.list()
        zeta = CyclotomicField(N).gen()
        iFT = [sum([S[i] * zeta**(-i * j) for i in J]) for j in J]
        if (J[0] not in ZZ) or F.base_ring().fraction_field() != QQ:
            raise NotImplementedError("Sorry this type of idft is not implemented yet.")
        return IndexedSequence(iFT, J) * (Integer(1) / N)

    def dct(self):
        """
        A discrete Cosine transform.

        EXAMPLES::

            sage: J = list(range(5))
            sage: A = [exp(-2*pi*i*I/5) for i in J]                                     # needs sage.symbolic
            sage: s = IndexedSequence(A, J)                                             # needs sage.symbolic
            sage: s.dct()                                                               # needs sage.symbolic
            Indexed sequence: [0, 1/16*(sqrt(5) + I*sqrt(-2*sqrt(5) + 10) + ...
            indexed by [0, 1, 2, 3, 4]
        """
        F = self.base_ring()      # elements must be coercible into RR
        try:
            pi = F.pi()
        except AttributeError:
            from sage.symbolic.constants import pi
            pi = F(pi)

        J = self.index_object()   # must be = range(N)
        N = len(J)
        S = self.list()
        PI = 2 * pi / N
        FT = [sum([S[i] * cos(PI * i * j) for i in J]) for j in J]
        return IndexedSequence(FT, J)

    def dst(self):
        """
        A discrete Sine transform.

        EXAMPLES::

            sage: J = list(range(5))
            sage: I = CC.0; pi = CC.pi()
            sage: A = [exp(-2*pi*i*I/5) for i in J]
            sage: s = IndexedSequence(A, J)

            sage: s.dst()        # discrete sine
            Indexed sequence: [0.000000000000000, 1.11022302462516e-16 - 2.50000000000000*I, ...]
            indexed by [0, 1, 2, 3, 4]
        """
        F = self.base_ring()      # elements must be coercible into RR
        try:
            pi = F.pi()
        except AttributeError:
            from sage.symbolic.constants import pi
            pi = F(pi)

        J = self.index_object()   # must be = range(N)
        N = len(J)
        S = self.list()
        PI = 2 * F(pi) / N
        FT = [sum([S[i] * sin(PI * i * j) for i in J]) for j in J]
        return IndexedSequence(FT, J)

    def convolution(self, other):
        r"""
        Convolves two sequences of the same length (automatically expands
        the shortest one by extending it by 0 if they have different lengths).

        If `\{a_n\}` and `\{b_n\}` are sequences indexed by `(n=0,1,...,N-1)`,
        extended by zero for all `n` in `\ZZ`, then the convolution is

        .. MATH::

             c_j = \sum_{i=0}^{N-1} a_i b_{j-i}.

        INPUT:

        - ``other`` -- a collection of elements of a ring with
          index set a finite abelian group (under `+`)

        OUTPUT: the Dirichlet convolution of ``self`` and ``other``

        EXAMPLES::

            sage: J = list(range(5))
            sage: A = [ZZ(1) for i in J]
            sage: B = [ZZ(1) for i in J]
            sage: s = IndexedSequence(A,J)
            sage: t = IndexedSequence(B,J)
            sage: s.convolution(t)
            [1, 2, 3, 4, 5, 4, 3, 2, 1]

        AUTHOR: David Joyner (2006-09)
        """
        S = self.list()
        T = other.list()
        I0 = self.index_object()
        J0 = other.index_object()
        F = self.base_ring()
        E = other.base_ring()
        if F != E:
            raise TypeError("IndexedSequences must have same base ring")
        if I0 != J0:
            raise TypeError("IndexedSequences must have same index set")
        M = len(S)
        N = len(T)
        if M < N:             # first, extend by 0 if necessary
            a = [S[i] for i in range(M)] + [F(0) for i in range(2 * N)]
            b = T + [E(0) for i in range(2 * M)]
        if M > N:             # python trick - a[-j] is really j from the *right*
            b = [T[i] for i in range(N)] + [E(0) for i in range(2 * M)]
            a = S + [F(0) for i in range(2 * M)]
        if M == N:              # so need only extend by 0 to the *right*
            a = S + [F(0) for i in range(2 * M)]
            b = T + [E(0) for i in range(2 * M)]
        N = max(M, N)
        return [sum([a[i] * b[j - i] for i in range(N)])
                for j in range(2 * N - 1)]

    def convolution_periodic(self, other):
        r"""
        Convolves two collections indexed by a ``range(...)`` of the same
        length (automatically expands the shortest one by extending it
        by 0 if they have different lengths).

        If `\{a_n\}` and `\{b_n\}` are sequences indexed by `(n=0,1,...,N-1)`,
        extended periodically for all `n` in `\ZZ`, then the convolution is

        .. MATH::

             c_j = \sum_{i=0}^{N-1} a_i b_{j-i}.

        INPUT:

        - ``other`` -- a sequence of elements of `\CC`, `\RR` or `\GF{q}`

        OUTPUT: the Dirichlet convolution of ``self`` and ``other``

        EXAMPLES::

            sage: I = list(range(5))
            sage: A = [ZZ(1) for i in I]
            sage: B = [ZZ(1) for i in I]
            sage: s = IndexedSequence(A,I)
            sage: t = IndexedSequence(B,I)
            sage: s.convolution_periodic(t)
            [5, 5, 5, 5, 5, 5, 5, 5, 5]

        AUTHOR: David Joyner (2006-09)
        """
        S = self.list()
        T = other.list()
        I = self.index_object()
        J = other.index_object()
        F = self.base_ring()
        E = other.base_ring()
        if F != E:
            raise TypeError("IndexedSequences must have same parent")
        if I != J:
            raise TypeError("IndexedSequences must have same index set")
        M = len(S)
        N = len(T)
        if M < N:           # first, extend by 0 if necessary
            a = [S[i] for i in range(M)] + [F(0) for i in range(N - M)]
            b = other
        if M > N:
            b = [T[i] for i in range(N)] + [E(0) for i in range(M - N)]
            a = self
        if M == N:
            a = S
            b = T
        N = max(M, N)
        return [sum([a[i] * b[(j - i) % N] for i in range(N)])
                for j in range(2 * N - 1)]

    def __mul__(self, other):
        """
        Implement scalar multiplication (on the right).

        EXAMPLES::

            sage: J = list(range(5))
            sage: A = [ZZ(1) for i in J]
            sage: s = IndexedSequence(A,J)
            sage: s.base_ring()
            Integer Ring
            sage: t = s*(1/3); t; t.base_ring()
            Indexed sequence: [1/3, 1/3, 1/3, 1/3, 1/3]
                indexed by [0, 1, 2, 3, 4]
            Rational Field
        """
        S = self.list()
        S1 = [S[i] * other for i in range(len(self.index_object()))]
        return IndexedSequence(S1, self.index_object())

    def __eq__(self, other):
        """
        Implement boolean equals.

        EXAMPLES::

            sage: J = list(range(5))
            sage: A = [ZZ(1) for i in J]
            sage: s = IndexedSequence(A,J)
            sage: t = s*(1/3)
            sage: t*3 == s
            1

        .. WARNING::

            ** elements are considered different if they differ
            by ``10^(-8)``, which is pretty arbitrary -- use with CAUTION!! **
        """
        if type(self) is not type(other):
            return False
        S = self.list()
        T = other.list()
        I = self.index_object()
        J = other.index_object()
        if I != J:
            return False
        for i in I:
            try:
                if abs(S[i] - T[i]) > 10**(-8):
                    # tests if they differ as reals  -- WHY 10^(-8)???
                    return False
            except TypeError:
                pass
        # if F != E:           # omitted this test since it
        #     return 0         # doesn't take into account coercions  -- WHY???
        return True

    def fft(self):
        """
        Wraps the gsl ``FastFourierTransform.forward()`` in
        :mod:`~sage.calculus.transforms.fft`.

        If the length is a power of 2 then this automatically uses the
        radix2 method. If the number of sample points in the input is
        a power of 2 then the wrapper for the GSL function
        ``gsl_fft_complex_radix2_forward()`` is automatically called.
        Otherwise, ``gsl_fft_complex_forward()`` is used.

        EXAMPLES::

            sage: J = list(range(5))
            sage: A = [RR(1) for i in J]
            sage: s = IndexedSequence(A,J)
            sage: t = s.fft(); t
            Indexed sequence: [5.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000]
                indexed by [0, 1, 2, 3, 4]
        """
        from sage.rings.cc import CC
        I = CC.gen()

        # elements must be coercible into RR
        J = self.index_object()   # must be = range(N)
        N = len(J)
        S = self.list()
        a = FastFourierTransform(N)
        for i in range(N):
            a[i] = S[i]
        a.forward_transform()
        return IndexedSequence([a[j][0] + I * a[j][1] for j in J], J)

    def ifft(self):
        """
        Implement the gsl ``FastFourierTransform.inverse`` in
        :mod:`~sage.calculus.transforms.fft`.

        If the number of sample points in the input is a power of 2
        then the wrapper for the GSL function
        ``gsl_fft_complex_radix2_inverse()`` is automatically called.
        Otherwise, ``gsl_fft_complex_inverse()`` is used.

        EXAMPLES::

            sage: J = list(range(5))
            sage: A = [RR(1) for i in J]
            sage: s = IndexedSequence(A,J)
            sage: t = s.fft(); t
            Indexed sequence: [5.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000]
                indexed by [0, 1, 2, 3, 4]
            sage: t.ifft()
            Indexed sequence: [1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000]
                indexed by [0, 1, 2, 3, 4]
            sage: t.ifft() == s
            1
        """
        from sage.rings.cc import CC
        I = CC.gen()

        # elements must be coercible into RR
        J = self.index_object()   # must be = range(N)
        N = len(J)
        S = self.list()
        a = FastFourierTransform(N)
        for i in range(N):
            a[i] = S[i]
        a.inverse_transform()
        return IndexedSequence([a[j][0] + I * a[j][1] for j in J], J)

    def dwt(self, other='haar', wavelet_k=2):
        r"""
        Wraps the gsl ``WaveletTransform.forward`` in :mod:`~sage.calculus.transforms.dwt`
        (written by Joshua Kantor). Assumes the length of the sample is a
        power of 2. Uses the GSL function ``gsl_wavelet_transform_forward()``.

        INPUT:

        - ``other`` -- the name of the type of wavelet; valid choices are:

          * ``'daubechies'``
          * ``'daubechies_centered'``
          * ``'haar'`` (default)
          * ``'haar_centered'``
          * ``'bspline'``
          * ``'bspline_centered'``

        - ``wavelet_k`` -- for daubechies wavelets, ``wavelet_k`` specifies a
          daubechie wavelet with `k/2` vanishing moments;
          `k = 4,6,...,20` for `k` even are the only ones implemented

          For Haar wavelets, ``wavelet_k`` must be 2.

          For bspline wavelets, ``wavelet_k`` equal to `103,105,202,204,
          206,208,301,305,307,309` will give biorthogonal B-spline wavelets
          of order `(i,j)` where ``wavelet_k`` equals `100 \cdot i + j`.

        The wavelet transform uses `J = \log_2(n)` levels.

        EXAMPLES::

            sage: J = list(range(8))
            sage: A = [RR(1) for i in J]
            sage: s = IndexedSequence(A,J)
            sage: t = s.dwt()
            sage: t        # slightly random output
            Indexed sequence: [2.82842712474999, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000]
                indexed by [0, 1, 2, 3, 4, 5, 6, 7]
        """
        from sage.rings.real_mpfr import RR
        # elements must be coercible into RR
        J = self.index_object()   # must be = range(N)
        N = len(J)             # must be 1 minus a power of 2
        S = self.list()
        if other == "haar" or other == "haar_centered":
            if wavelet_k in [2]:
                a = WaveletTransform(N, other, wavelet_k)
            else:
                raise ValueError("wavelet_k must be = 2")
        if other == "daubechies" or other == "daubechies_centered":
            if wavelet_k in [4, 6, 8, 10, 12, 14, 16, 18, 20]:
                a = WaveletTransform(N, other, wavelet_k)
            else:
                raise ValueError("wavelet_k must be in {4,6,8,10,12,14,16,18,20}")
        if other == "bspline" or other == "bspline_centered":
            if wavelet_k in [103, 105, 202, 204, 206, 208, 301, 305, 307, 309]:
                a = WaveletTransform(N, other, 103)
            else:
                raise ValueError("wavelet_k must be in {103,105,202,204,206,208,301,305,307,309}")
        for i in range(N):
            a[i] = S[i]
        a.forward_transform()
        return IndexedSequence([RR(a[j]) for j in J], J)

    def idwt(self, other='haar', wavelet_k=2):
        r"""
        Implement the gsl ``WaveletTransform.backward()`` in
        :mod:`~sage.calculus.transforms.dwt`.

        Assumes the length of the sample is a power of 2. Uses the
        GSL function ``gsl_wavelet_transform_backward()``.

        INPUT:

        - ``other`` -- must be one of the following:

          * ``'haar'``
          * ``'daubechies'``
          * ``'daubechies_centered'``
          * ``'haar_centered'``
          * ``'bspline'``
          * ``'bspline_centered'``

        - ``wavelet_k`` -- for daubechies wavelets, ``wavelet_k`` specifies a
          daubechie wavelet with `k/2` vanishing moments;
          `k = 4,6,...,20` for `k` even are the only ones implemented

          For Haar wavelets, ``wavelet_k`` must be 2.

          For bspline wavelets, ``wavelet_k`` equal to `103,105,202,204,
          206,208,301,305,307,309` will give biorthogonal B-spline wavelets
          of order `(i,j)` where ``wavelet_k`` equals `100 \cdot i + j`.

        EXAMPLES::

            sage: J = list(range(8))
            sage: A = [RR(1) for i in J]
            sage: s = IndexedSequence(A,J)
            sage: t = s.dwt()
            sage: t            # random arch dependent output
            Indexed sequence: [2.82842712474999, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000]
                indexed by [0, 1, 2, 3, 4, 5, 6, 7]
            sage: t.idwt()                  # random arch dependent output
            Indexed sequence: [1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000]
                indexed by [0, 1, 2, 3, 4, 5, 6, 7]
            sage: t.idwt() == s
            True
            sage: J = list(range(16))
            sage: A = [RR(1) for i in J]
            sage: s = IndexedSequence(A,J)
            sage: t = s.dwt("bspline", 103)
            sage: t   # random arch dependent output
            Indexed sequence: [4.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000]
                indexed by [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
            sage: t.idwt("bspline", 103) == s
            True
        """
        from sage.rings.real_mpfr import RR
        # elements must be coercible into RR
        J = self.index_object()   # must be = range(N)
        N = len(J)             # must be 1 minus a power of 2
        S = self.list()
        k = wavelet_k
        if other == "haar" or other == "haar_centered":
            if k in [2]:
                a = WaveletTransform(N, other, wavelet_k)
            else:
                raise ValueError("wavelet_k must be = 2")
        if other == "daubechies" or other == "daubechies_centered":
            if k in [4, 6, 8, 10, 12, 14, 16, 18, 20]:
                a = WaveletTransform(N, other, wavelet_k)
            else:
                raise ValueError("wavelet_k must be in {4,6,8,10,12,14,16,18,20}")
        if other == "bspline" or other == "bspline_centered":
            if k in [103, 105, 202, 204, 206, 208, 301, 305, 307, 309]:
                a = WaveletTransform(N, other, 103)
            else:
                raise ValueError("wavelet_k must be in {103,105,202,204,206,208,301,305,307,309}")
        for i in range(N):
            a[i] = S[i]
        a.backward_transform()
        return IndexedSequence([RR(a[j]) for j in J], J)
