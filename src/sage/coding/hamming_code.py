# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Hamming codes

Given an integer `r` and a field `F`, such that `F=\GF{q}`, the `[n, k, d]` code
with length `n=\frac{q^{r}-1}{q-1}`, dimension `k=\frac{q^{r}-1}{q-1} - r` and
minimum distance `d=3` is called the Hamming Code of order `r`.

REFERENCES:

- [Rot2006]_
"""
# ****************************************************************************
#       Copyright (C) 2016 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .linear_code import AbstractLinearCode
from sage.matrix.matrix_space import MatrixSpace
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer


class HammingCode(AbstractLinearCode):
    r"""
    Representation of a Hamming code.

    INPUT:

    - ``base_field`` -- the base field over which ``self`` is defined

    - ``order`` -- the order of ``self``

    EXAMPLES::

        sage: C = codes.HammingCode(GF(7), 3)
        sage: C
        [57, 54] Hamming Code over GF(7)
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, order):
        r"""
        TESTS:

        If ``base_field`` is not a finite field, an exception is raised::

            sage: codes.HammingCode(RR, 3)
            Traceback (most recent call last):
            ...
            ValueError: base_field has to be a finite field

        If ``order`` is not a Sage Integer or a Python int, an exception is raised::

            sage: codes.HammingCode(GF(3), 3.14)
            Traceback (most recent call last):
            ...
            ValueError: order has to be a Sage Integer or a Python int
        """
        if not base_field.is_finite():
            raise ValueError("base_field has to be a finite field")
        if not isinstance(order, (Integer, int)):
            raise ValueError("order has to be a Sage Integer or a Python int")

        q = base_field.order()
        length = Integer((q ** order - 1) / (q - 1))
        super().__init__(base_field, length, "Systematic", "Syndrome")
        self._dimension = length - order

    def __eq__(self, other):
        r"""
        Test equality of Hamming Code objects.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(7), 3)
            sage: C2 = codes.HammingCode(GF(7), 3)
            sage: C1 == C2
            True
        """
        return isinstance(other, HammingCode)\
                and self.length() == other.length()\
                and self.dimension() == other.dimension()

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(7), 3)
            sage: C2 = codes.HammingCode(GF(7), 3)
            sage: hash(C1) == hash(C2)
            True
        """
        return hash((self.length(), self.dimension()))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(7), 3)
            sage: C
            [57, 54] Hamming Code over GF(7)
        """
        return "[%s, %s] Hamming Code over GF(%s)"\
                % (self.length(), self.dimension(), self.base_field().cardinality())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(7), 3)
            sage: latex(C)
            [57, 54] \textnormal{ Hamming Code over }\Bold{F}_{7}
        """
        return "[%s, %s] \\textnormal{ Hamming Code over }%s"\
                % (self.length(), self.dimension(), self.base_field()._latex_())

    @cached_method
    def parity_check_matrix(self):
        r"""
        Return a parity check matrix of ``self``.

        The construction of the parity check matrix in case ``self``
        is not a binary code is not really well documented.
        Regarding the choice of projective geometry, one might check:

        - the note over section 2.3 in [Rot2006]_, pages 47-48
        - the dedicated paragraph in [HP2003]_, page 30

        EXAMPLES::

            sage: C = codes.HammingCode(GF(3), 3)
            sage: C.parity_check_matrix()
            [1 0 1 1 0 1 0 1 1 1 0 1 1]
            [0 1 1 2 0 0 1 1 2 0 1 1 2]
            [0 0 0 0 1 1 1 1 1 2 2 2 2]
        """
        n = self.length()
        F = self.base_field()
        m = n - self.dimension()
        MS = MatrixSpace(F, n, m)
        X = ProjectiveSpace(m - 1, F)
        PFn = [list(p) for p in X.point_set(F).points()]

        H = MS(PFn).transpose()
        H = H[::-1, :]
        H.set_immutable()
        return H

    def minimum_distance(self):
        r"""
        Return the minimum distance of ``self``.

        It is always 3 as ``self`` is a Hamming Code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(7), 3)
            sage: C.minimum_distance()
            3
        """
        return 3
