# sage_setup: distribution = sagemath-modules
r"""
Orthogonal Linear Groups

The general orthogonal group `GO(n,R)` consists of all `n \times n`
matrices over the ring `R` preserving an `n`-ary positive definite
quadratic form. In cases where there are multiple non-isomorphic
quadratic forms, additional data needs to be specified to
disambiguate. The special orthogonal group is the normal subgroup of
matrices of determinant one.

In characteristics different from 2, a quadratic form is equivalent to
a bilinear symmetric form. Furthermore, over the real numbers a
positive definite quadratic form is equivalent to the diagonal
quadratic form, equivalent to the bilinear symmetric form defined by
the identity matrix. Hence, the orthogonal group `GO(n,\RR)` is the
group of orthogonal matrices in the usual sense.

In the case of a finite field and if the degree `n` is even, then
there are two inequivalent quadratic forms and a third parameter ``e``
must be specified to disambiguate these two possibilities. The index
of `SO(e,d,q)` in `GO(e,d,q)` is `2` if `q` is odd, but `SO(e,d,q) =
GO(e,d,q)` if `q` is even.)

.. WARNING::

   GAP and Sage use different notations:

   * GAP notation: The optional ``e`` comes first, that is, ``GO([e,]
     d, q)``, ``SO([e,] d, q)``.

   * Sage notation: The optional ``e`` comes last, the standard Python
     convention: ``GO(d, GF(q), e=0)``, ``SO(d, GF(q), e=0)``.

EXAMPLES::

    sage: GO(3,7)
    General Orthogonal Group of degree 3 over Finite Field of size 7

    sage: G = SO(4, GF(7), 1); G
    Special Orthogonal Group of degree 4 and form parameter 1
     over Finite Field of size 7

    sage: # needs sage.libs.gap
    sage: G.random_element()   # random
    [4 3 5 2]
    [6 6 4 0]
    [0 4 6 0]
    [4 4 5 1]

TESTS::

    sage: G = GO(3, GF(5))
    sage: latex(G)
    \text{GO}_{3}(\Bold{F}_{5})
    sage: G = SO(3, GF(5))
    sage: latex(G)
    \text{SO}_{3}(\Bold{F}_{5})
    sage: G = SO(4, GF(5), 1)
    sage: latex(G)
    \text{SO}_{4}(\Bold{F}_{5}, +)

AUTHORS:

- David Joyner (2006-03): initial version

- David Joyner (2006-05): added examples, _latex_, __str__, gens,
  as_matrix_group

- William Stein (2006-12-09): rewrite

- Volker Braun (2013-1) port to new Parent, libGAP, extreme refactoring.

- Sebastian Oehms (2018-8) add
  :meth:`~sage.groups.matrix_gps.orthogonal.OrthogonalMatrixGroup_generic.invariant_form`
  (as alias), ``_OG``, option for user defined invariant bilinear form,
  and bug-fix in cmd-string for calling GAP (see :issue:`26028`)
"""

# ****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#                     2009 Mike Hansen
#                     2013 Volker Braun <vbraun.name@gmail.com>
#                     2018 Sebastian Oehms
#                     2018 Travis Scrimshaw
#                     2023 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.misc.latex import latex
from sage.misc.cachefunc import cached_method
from sage.groups.matrix_gps.named_group import (
    normalize_args_vectorspace, normalize_args_invariant_form,
    NamedMatrixGroup_generic)


def normalize_args_e(degree, ring, e):
    """
    Normalize the arguments that relate the choice of quadratic form
    for special orthogonal groups over finite fields.

    INPUT:

    - ``degree`` -- integer; the degree of the affine group, that is,
      the dimension of the affine space the group is acting on

    - ``ring`` -- a ring; the base ring of the affine space

    - ``e`` -- integer; one of `+1`, `0`, `-1`.  Only relevant for
      finite fields and if the degree is even. A parameter that
      distinguishes inequivalent invariant forms.

    OUTPUT: the integer ``e`` with values required by GAP

    TESTS::

        sage: from sage.groups.matrix_gps.orthogonal import normalize_args_e
        sage: normalize_args_e(2, GF(3), +1)
        1
        sage: normalize_args_e(3, GF(3), 0)
        0
        sage: normalize_args_e(3, GF(3), +1)
        0
        sage: normalize_args_e(2, GF(3), 0)
        Traceback (most recent call last):
        ...
        ValueError: must have e=-1 or e=1 for even degree
    """
    if isinstance(ring, FiniteField) and degree % 2 == 0:
        if e not in (-1, +1):
            raise ValueError('must have e=-1 or e=1 for even degree')
    else:
        e = 0
    return ZZ(e)


###############################################################################
# Orthogonal Group: common Code for both GO and SO
###############################################################################

def _OG(n, R, special, e=0, var='a', invariant_form=None):
    r"""
    This function is commonly used by the functions GO and SO to avoid
    unnecessarily duplicated code. For documentation and examples see
    the individual functions.

    TESTS:

    Check that :issue:`26028` is fixed::

        sage: GO(3,25).order()  # indirect doctest                                      # needs sage.libs.gap sage.rings.finite_rings
        31200

    Check that :issue:`28054` is fixed::

        sage: G = SO(2, GF(3), -1)
        sage: m = G.invariant_form()
        sage: G2 = SO(2, GF(3), 1, invariant_form=m)
        Traceback (most recent call last):
        ...
        NotImplementedError: invariant_form for finite groups is fixed by GAP
    """
    prefix = 'General'
    ltx_prefix = 'G'
    if special:
        prefix = 'Special'
        ltx_prefix = 'S'

    degree, ring = normalize_args_vectorspace(n, R, var=var)
    e = normalize_args_e(degree, ring, e)

    if invariant_form is not None:
        if isinstance(ring, FiniteField):
            raise NotImplementedError("invariant_form for finite groups is fixed by GAP")

    if e == 0:
        if invariant_form is not None:
            invariant_form = normalize_args_invariant_form(ring, degree, invariant_form)
            if not invariant_form.is_symmetric():
                raise ValueError("invariant_form must be symmetric")

            try:
                if invariant_form.is_positive_definite():
                    inserted_text = "with respect to positive definite symmetric form"
                else:
                    inserted_text = "with respect to non positive definite symmetric form"
            except ValueError:
                inserted_text = "with respect to symmetric form"

            name = '{0} Orthogonal Group of degree {1} over {2} {3}\n{4}'.format(
                prefix, degree, ring, inserted_text, invariant_form)
            ltx = r'\text{{{0}O}}_{{{1}}}({2})\text{{ {3} }}{4}'.format(
                ltx_prefix, degree, latex(ring), inserted_text,
                latex(invariant_form))
        else:
            name = '{0} Orthogonal Group of degree {1} over {2}'.format(prefix, degree, ring)
            ltx = r'\text{{{0}O}}_{{{1}}}({2})'.format(ltx_prefix, degree, latex(ring))
    else:
        name = '{0} Orthogonal Group of degree {1} and form parameter {2} over {3}'.format(prefix, degree, e, ring)
        ltx = r'\text{{{0}O}}_{{{1}}}({2}, {3})'.format(ltx_prefix, degree,
                                                        latex(ring),
                                                        '+' if e == 1 else '-')

    if isinstance(ring, FiniteField):
        try:
            from .orthogonal_gap import OrthogonalMatrixGroup_gap
        except ImportError:
            pass
        else:
            cmd = '{0}O({1}, {2}, {3})'.format(ltx_prefix, e, degree, ring.order())
            return OrthogonalMatrixGroup_gap(degree, ring, False, name, ltx, cmd)

    return OrthogonalMatrixGroup_generic(degree, ring, False, name, ltx, invariant_form=invariant_form)


########################################################################
# General Orthogonal Group
########################################################################

def GO(n, R, e=0, var='a', invariant_form=None):
    r"""
    Return the general orthogonal group.

    The general orthogonal group `GO(n,R)` consists of all `n \times n`
    matrices over the ring `R` preserving an `n`-ary positive definite
    quadratic form. In cases where there are multiple non-isomorphic
    quadratic forms, additional data needs to be specified to
    disambiguate.

    In the case of a finite field and if the degree `n` is even, then
    there are two inequivalent quadratic forms and a third parameter
    ``e`` must be specified to disambiguate these two possibilities.

    .. NOTE::

        This group is also available via ``groups.matrix.GO()``.

    INPUT:

    - ``n`` -- integer; the degree

    - ``R`` -- ring or an integer; if an integer is specified, the
      corresponding finite field is used

    - ``e`` -- ``+1`` or ``-1``, and ignored by default; only relevant
      for finite fields and if the degree is even: a parameter that
      distinguishes inequivalent invariant forms

    - ``var`` -- (default: ``'a'``) variable used to
      represent generator of the finite field, if needed

    - ``invariant_form`` -- (optional) instances being accepted by
      the matrix-constructor which define a `n \times n` square matrix
      over ``R`` describing the symmetric form to be kept invariant
      by the orthogonal group; the form is checked to be
      non-degenerate and symmetric but not to be positive definite

    OUTPUT: the general orthogonal group of given degree, base ring, and
    choice of invariant form

    EXAMPLES::

        sage: GO(3, GF(7))
        General Orthogonal Group of degree 3 over Finite Field of size 7

        sage: # needs sage.libs.gap
        sage: GO(3, GF(7)).order()
        672
        sage: GO(3, GF(7)).gens()
        (
        [3 0 0]  [0 1 0]
        [0 5 0]  [1 6 6]
        [0 0 1], [0 2 1]
        )

    Using the ``invariant_form`` option::

        sage: m = matrix(QQ, 3, 3, [[0, 1, 0], [1, 0, 0], [0, 0, 3]])
        sage: GO3  = GO(3, QQ)
        sage: GO3m = GO(3, QQ, invariant_form=m)
        sage: GO3 == GO3m
        False
        sage: GO3.invariant_form()
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: GO3m.invariant_form()
        [0 1 0]
        [1 0 0]
        [0 0 3]
        sage: pm = Permutation([2,3,1]).to_matrix()
        sage: g = GO3(pm); g in GO3; g
        True
        [0 0 1]
        [1 0 0]
        [0 1 0]
        sage: GO3m(pm)
        Traceback (most recent call last):
        ...
        TypeError: matrix must be orthogonal with respect to the symmetric form
        [0 1 0]
        [1 0 0]
        [0 0 3]

        sage: GO(3,3, invariant_form=[[1,0,0], [0,2,0], [0,0,1]])
        Traceback (most recent call last):
        ...
        NotImplementedError: invariant_form for finite groups is fixed by GAP
        sage: 5 + 5
        10
        sage: R.<x> = ZZ[]
        sage: GO(2, R, invariant_form=[[x,0], [0,1]])
        General Orthogonal Group of degree 2 over
         Univariate Polynomial Ring in x over Integer Ring with respect to symmetric form
        [x 0]
        [0 1]

    TESTS::

        sage: # needs sage.libs.gap
        sage: TestSuite(GO3).run()
        sage: groups.matrix.GO(2, 3, e=-1)                                              # needs sage.rings.finite_rings
        General Orthogonal Group of degree 2 and form parameter -1 over Finite Field of size 3
    """
    return _OG(n, R, False, e=e, var=var, invariant_form=invariant_form)


########################################################################
# Special Orthogonal Group
########################################################################

def SO(n, R, e=None, var='a', invariant_form=None):
    r"""
    Return the special orthogonal group.

    The special orthogonal group `GO(n,R)` consists of all `n \times n`
    matrices with determinant one over the ring `R` preserving an
    `n`-ary positive definite quadratic form. In cases where there are
    multiple non-isomorphic quadratic forms, additional data needs to
    be specified to disambiguate.

    .. NOTE::

        This group is also available via ``groups.matrix.SO()``.

    INPUT:

    - ``n`` -- integer; the degree

    - ``R`` -- ring or an integer; if an integer is specified, the
      corresponding finite field is used

    - ``e`` -- ``+1`` or ``-1``, and ignored by default; only relevant
      for finite fields and if the degree is even: a parameter that
      distinguishes inequivalent invariant forms

    - ``var`` -- (default: ``'a'``) variable used to
      represent generator of the finite field, if needed

    - ``invariant_form`` -- (optional) instances being accepted by
      the matrix-constructor which define a `n \times n` square matrix
      over ``R`` describing the symmetric form to be kept invariant
      by the orthogonal group; the form is checked to be
      non-degenerate and symmetric but not to be positive definite

    OUTPUT: the special orthogonal group of given degree, base ring, and choice
    of invariant form

    EXAMPLES::

        sage: G = SO(3,GF(5)); G
        Special Orthogonal Group of degree 3 over Finite Field of size 5

        sage: # needs sage.libs.gap
        sage: G = SO(3,GF(5))
        sage: G.gens()
        (
        [2 0 0]  [3 2 3]  [1 4 4]
        [0 3 0]  [0 2 0]  [4 0 0]
        [0 0 1], [0 3 1], [2 0 4]
        )
        sage: G.as_matrix_group()
        Matrix group over Finite Field of size 5 with 3 generators (
        [2 0 0]  [3 2 3]  [1 4 4]
        [0 3 0]  [0 2 0]  [4 0 0]
        [0 0 1], [0 3 1], [2 0 4]
        )

    Using the ``invariant_form`` option::

        sage: # needs sage.rings.number_field
        sage: CF3 = CyclotomicField(3); e3 = CF3.gen()
        sage: m = matrix(CF3, 3, 3, [[1,e3,0], [e3,2,0], [0,0,1]])
        sage: SO3  = SO(3, CF3)
        sage: SO3m = SO(3, CF3, invariant_form=m)
        sage: SO3 == SO3m
        False
        sage: SO3.invariant_form()
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: SO3m.invariant_form()
        [    1 zeta3     0]
        [zeta3     2     0]
        [    0     0     1]
        sage: pm = Permutation([2,3,1]).to_matrix()
        sage: g = SO3(pm); g in SO3; g
        True
        [0 0 1]
        [1 0 0]
        [0 1 0]
        sage: SO3m(pm)
        Traceback (most recent call last):
        ...
        TypeError: matrix must be orthogonal with respect to the symmetric form
        [    1 zeta3     0]
        [zeta3     2     0]
        [    0     0     1]

        sage: SO(3, 5, invariant_form=[[1,0,0], [0,2,0], [0,0,3]])
        Traceback (most recent call last):
        ...
        NotImplementedError: invariant_form for finite groups is fixed by GAP
        sage: 5+5
        10

    TESTS::

        sage: TestSuite(SO3m).run()                                                     # needs sage.rings.number_field
        sage: groups.matrix.SO(2, 3, e=1)
        Special Orthogonal Group of degree 2 and form parameter 1 over Finite Field of size 3
    """
    return _OG(n, R, True, e=e, var=var, invariant_form=invariant_form)


########################################################################
# Orthogonal Group class
########################################################################

class OrthogonalMatrixGroup_generic(NamedMatrixGroup_generic):
    r"""
    General Orthogonal Group over arbitrary rings.

    EXAMPLES::

        sage: G = GO(3, GF(7)); G
        General Orthogonal Group of degree 3 over Finite Field of size 7
        sage: latex(G)
        \text{GO}_{3}(\Bold{F}_{7})

        sage: G = SO(3, GF(5));  G
        Special Orthogonal Group of degree 3 over Finite Field of size 5
        sage: latex(G)
        \text{SO}_{3}(\Bold{F}_{5})

        sage: # needs sage.rings.number_field
        sage: CF3 = CyclotomicField(3); e3 = CF3.gen()
        sage: m = matrix(CF3, 3,3, [[1,e3,0],[e3,2,0],[0,0,1]])
        sage: G = SO(3, CF3, invariant_form=m)
        sage: latex(G)
        \text{SO}_{3}(\Bold{Q}(\zeta_{3}))\text{ with respect to non positive definite symmetric form }\left(\begin{array}{rrr}
        1 & \zeta_{3} & 0 \\
        \zeta_{3} & 2 & 0 \\
        0 & 0 & 1
        \end{array}\right)
    """

    @cached_method
    def invariant_bilinear_form(self):
        """
        Return the symmetric bilinear form preserved by ``self``.

        OUTPUT: a matrix

        EXAMPLES::

            sage: # needs sage.libs.gap
            sage: GO(2,3,+1).invariant_bilinear_form()
            [0 1]
            [1 0]
            sage: GO(2,3,-1).invariant_bilinear_form()
            [2 1]
            [1 1]
            sage: G = GO(4, QQ)
            sage: G.invariant_bilinear_form()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: GO3m = GO(3, QQ, invariant_form=(1,0,0, 0,2,0, 0,0,3))
            sage: GO3m.invariant_bilinear_form()
            [1 0 0]
            [0 2 0]
            [0 0 3]

        TESTS::

            sage: GO3m.invariant_form()                                                 # needs sage.libs.gap
            [1 0 0]
            [0 2 0]
            [0 0 3]
        """
        if self._invariant_form is not None:
            return self._invariant_form

        from sage.matrix.constructor import identity_matrix
        m = identity_matrix(self.base_ring(), self.degree())
        m.set_immutable()
        return m

    invariant_quadratic_form = invariant_bilinear_form
    # this is identical in the generic case

    invariant_form = invariant_bilinear_form
    # alias (analogues to symplectic and unitary cases)

    def _check_matrix(self, x, *args):
        """
        Check whether the matrix ``x`` is orthogonal.

        See :meth:`~sage.groups.matrix_gps.matrix_group._check_matrix`
        for details.

        EXAMPLES::

            sage: G = GO(4, GF(5), +1)
            sage: G._check_matrix(G.an_element().matrix())
        """
        if self._special and x.determinant() != 1:
            raise TypeError('matrix must have determinant one')
        F = self.invariant_bilinear_form()
        if x * F * x.transpose() != F:
            if F == self.one().matrix():
                raise TypeError('matrix must be orthogonal')
            else:
                raise TypeError('matrix must be orthogonal with respect to the symmetric form\n%s' % (F))
        # TODO: check that quadratic form is preserved in characteristic two
