# sage.doctest: needs sage.libs.pari
"""
Module of supersingular points

The module of divisors on the modular curve `X_0(N)` over `F_p` supported at supersingular points.

EXAMPLES::

    sage: x = SupersingularModule(389)
    sage: m = x.T(2).matrix()
    sage: a = m.change_ring(GF(97))
    sage: D = a.decomposition()
    sage: D[:3]
    [(Vector space of degree 33 and dimension 1 over Finite Field of size 97
      Basis matrix:
      [ 0  0  0  1 96 96  1  0 95  1  1  1  1 95  2 96  0  0 96  0 96  0 96  2 96 96  0  1  0  2  1 95  0],
      True),
     (Vector space of degree 33 and dimension 1 over Finite Field of size 97
      Basis matrix:
      [ 0  1 96 16 75 22 81  0  0 17 17 80 80  0  0 74 40  1 16 57 23 96 81  0 74 23  0 24  0  0 73  0  0],
      True),
     (Vector space of degree 33 and dimension 1 over Finite Field of size 97
      Basis matrix:
      [ 0  1 96 90 90  7  7  0  0 91  6  6 91  0  0 91  0 13  7  0  6 84 90  0  6 91  0 90  0  0  7  0  0],
      True)]
    sage: len(D)
    9

We compute a Hecke operator on a space of huge dimension!::

    sage: X = SupersingularModule(next_prime(10000))
    sage: t = X.T(2).matrix()            # long time (21s on sage.math, 2011)
    sage: t.nrows()                      # long time
    835

TESTS::

    sage: X = SupersingularModule(389)
    sage: T = X.T(2).matrix().change_ring(QQ)
    sage: d = T.decomposition()
    sage: len(d)
    6
    sage: [a[0].dimension() for a in d]
    [1, 1, 2, 3, 6, 20]
    sage: loads(dumps(X)) == X
    True
    sage: loads(dumps(d)) == d
    True

AUTHORS:

- William Stein

- David Kohel

- Iftikhar Burhanuddin
"""

# ****************************************************************************
#       Copyright (C) 2004, 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 Iftikhar Burhanuddin <burhanud@usc.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import kronecker, next_prime
from sage.categories.fields import Fields
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.lazy_import import lazy_import
from sage.modular.arithgroup.all import Gamma0
from sage.modular.hecke.module import HeckeModule_free_module
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.richcmp import richcmp_method, richcmp

lazy_import('sage.libs.pari', 'pari')


ZZy = PolynomialRing(ZZ, 'y')


def Phi2_quad(J3, ssJ1, ssJ2):
    r"""
    Return a certain quadratic polynomial over a finite
    field in indeterminate J3.

    The roots of the polynomial along with ssJ1 are the
    neighboring/2-isogenous supersingular j-invariants of ssJ2.

    INPUT:

    - ``J3`` -- indeterminate of a univariate polynomial ring defined over a finite
      field with p^2 elements where p is a prime number

    - ``ssJ2``, ``ssJ2`` -- supersingular j-invariants over the finite field

    OUTPUT: polynomial; defined over the finite field

    EXAMPLES:

    The following code snippet produces a factor of the modular polynomial
    `\Phi_{2}(x,j_{in})`, where `j_{in}` is a supersingular j-invariant
    defined over the finite field with `37^2` elements::

        sage: F = GF(37^2, 'a')
        sage: X = PolynomialRing(F, 'x').gen()
        sage: j_in = supersingular_j(F)
        sage: poly = sage.modular.ssmod.ssmod.Phi_polys(2,X,j_in)
        sage: poly.roots()
        [(8, 1), (27*a + 23, 1), (10*a + 20, 1)]
        sage: sage.modular.ssmod.ssmod.Phi2_quad(X, F(8), j_in)
        x^2 + 31*x + 31

    .. NOTE::

        Given a root (j1,j2) to the polynomial `Phi_2(J1,J2)`, the pairs
        (j2,j3) not equal to (j2,j1) which solve `Phi_2(j2,j3)` are roots of
        the quadratic equation:

        J3^2 + (-j2^2 + 1488*j2 + (j1 - 162000))*J3 + (-j1 + 1488)*j2^2 +
        (1488*j1 + 40773375)*j2 + j1^2 - 162000*j1 + 8748000000

        This will be of use to extend the 2-isogeny graph, once the initial
        three roots are determined for `Phi_2(J1,J2)`.

    AUTHORS:

    - David Kohel -- kohel@maths.usyd.edu.au

    - Iftikhar Burhanuddin -- burhanud@usc.edu
    """
    ssJ1_pow2 = ssJ1**2
    ssJ2_pow2 = ssJ2**2

    return J3.parent()([(-ssJ1 + 1488) * ssJ2_pow2
                        + (1488 * ssJ1 + 40773375) * ssJ2
                        + ssJ1_pow2 - 162000 * ssJ1 + 8748000000,
                        -ssJ2_pow2 + 1488 * ssJ2 + (ssJ1 - 162000),
                        1])


def Phi_polys(L, x, j):
    r"""
    Return a certain polynomial of degree `L+1` in the
    indeterminate x over a finite field.

    The roots of the **modular** polynomial `\Phi(L, x, j)` are the
    `L`-isogenous supersingular j-invariants of j.

    INPUT:

    - ``L`` -- integer

    - ``x`` -- indeterminate of a univariate polynomial ring defined over a
      finite field with p^2 elements, where p is a prime number

    - ``j`` -- supersingular j-invariant over the finite field

    OUTPUT: polynomial; defined over the finite field

    EXAMPLES:

    The following code snippet produces the modular polynomial
    `\Phi_{L}(x,j_{in})`, where `j_{in}` is a supersingular j-invariant
    defined over the finite field with `7^2` elements::

        sage: F = GF(7^2, 'a')
        sage: X = PolynomialRing(F, 'x').gen()
        sage: j_in = supersingular_j(F)
        sage: sage.modular.ssmod.ssmod.Phi_polys(2,X,j_in)
        x^3 + 3*x^2 + 3*x + 1
        sage: sage.modular.ssmod.ssmod.Phi_polys(3,X,j_in)
        x^4 + 4*x^3 + 6*x^2 + 4*x + 1
        sage: sage.modular.ssmod.ssmod.Phi_polys(5,X,j_in)
        x^6 + 6*x^5 + x^4 + 6*x^3 + x^2 + 6*x + 1
        sage: sage.modular.ssmod.ssmod.Phi_polys(7,X,j_in)
        x^8 + x^7 + x + 1
        sage: sage.modular.ssmod.ssmod.Phi_polys(11,X,j_in)
        x^12 + 5*x^11 + 3*x^10 + 3*x^9 + 5*x^8 + x^7 + x^5 + 5*x^4 + 3*x^3 + 3*x^2 + 5*x + 1
        sage: sage.modular.ssmod.ssmod.Phi_polys(13,X,j_in)
        x^14 + 2*x^7 + 1
    """
    r = 0
    for pol in pari.polmodular(L).Vec():
        r = r * x + ZZy(pol)(j)
    return r


def dimension_supersingular_module(prime, level=1):
    r"""
    Return the dimension of the Supersingular module, which is
    equal to the dimension of the space of modular forms of weight `2`
    and conductor equal to ``prime`` times ``level``.

    INPUT:

    - ``prime`` -- integer; prime

    - ``level`` -- integer; positive

    OUTPUT: dimension; integer, nonnegative

    EXAMPLES:

    The code below computes the dimensions of Supersingular modules
    with level=1 and prime = 7, 15073 and 83401::

        sage: dimension_supersingular_module(7)
        1

        sage: dimension_supersingular_module(15073)
        1256

        sage: dimension_supersingular_module(83401)
        6950

    .. NOTE::

        The case of level > 1 has not been implemented yet.

    AUTHORS:

    - David Kohel -- kohel@maths.usyd.edu.au

    - Iftikhar Burhanuddin - burhanud@usc.edu
    """
    if not Integer(prime).is_prime():
        raise ValueError(f"{prime} is not a prime")

    if level == 1:
        return Gamma0(prime).dimension_modular_forms(2)

    # list of genus(X_0(level)) equal to zero
    # elif (level in [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 18, 25]):
    # compute basis

    raise NotImplementedError


def supersingular_D(prime):
    r"""
    Return a fundamental discriminant `D` of an
    imaginary quadratic field, where the given prime does not split.

    See Silverman's Advanced Topics in the Arithmetic of Elliptic
    Curves, page 184, exercise 2.30(d).

    INPUT:

    - ``prime`` -- integer, prime

    OUTPUT: d; integer, negative

    EXAMPLES:

    These examples return *supersingular discriminants* for 7,
    15073 and 83401::

        sage: supersingular_D(7)
        -4

        sage: supersingular_D(15073)
        -15

        sage: supersingular_D(83401)
        -7

    AUTHORS:

    - David Kohel - kohel@maths.usyd.edu.au

    - Iftikhar Burhanuddin - burhanud@usc.edu
    """
    if not Integer(prime).is_prime():
        raise ValueError("%s is not a prime" % prime)

    # Making picking D more intelligent
    D = -1
    while True:
        Dmod4 = D % 4
        if Dmod4 in (0, 1) and kronecker(D, prime) != 1:
            return D
        D -= 1


def supersingular_j(FF):
    r"""
    Return a supersingular j-invariant over the finite
    field FF.

    INPUT:

    - ``FF`` -- finite field with p^2 elements, where p is a prime number

    OUTPUT:

    - finite field element -- a supersingular j-invariant
      defined over the finite field FF

    EXAMPLES:

    The following examples calculate supersingular j-invariants for a
    few finite fields::

        sage: supersingular_j(GF(7^2, 'a'))
        6

    Observe that in this example the j-invariant is not defined over
    the prime field::

        sage: supersingular_j(GF(15073^2, 'a'))
        4443*a + 13964
        sage: supersingular_j(GF(83401^2, 'a'))
        67977

    AUTHORS:

    - David Kohel -- kohel@maths.usyd.edu.au

    - Iftikhar Burhanuddin -- burhanud@usc.edu
    """
    if FF not in Fields().Finite():
        raise ValueError("%s is not a finite field" % FF)
    prime = FF.characteristic()
    if not Integer(prime).is_prime():
        raise ValueError("%s is not a prime" % prime)
    if FF.cardinality() != Integer(prime**2):
        raise ValueError("%s is not a quadratic extension" % FF)
    if kronecker(-1, prime) != 1:
        j_invss = 1728                 # (2^2 * 3)^3
    elif kronecker(-2, prime) != 1:
        j_invss = 8000                 # (2^2 * 5)^3
    elif kronecker(-3, prime) != 1:
        j_invss = 0                    # 0^3
    elif kronecker(-7, prime) != 1:
        j_invss = 16581375             # (3 * 5 * 17)^3
    elif kronecker(-11, prime) != 1:
        j_invss = -32768               # -(2^5)^3
    elif kronecker(-19, prime) != 1:
        j_invss = -884736              # -(2^5 * 3)^3
    elif kronecker(-43, prime) != 1:
        j_invss = -884736000           # -(2^6 * 3 * 5)^3
    elif kronecker(-67, prime) != 1:
        j_invss = -147197952000        # -(2^5 * 3 * 5 * 11)^3
    elif kronecker(-163, prime) != 1:
        j_invss = -262537412640768000  # -(2^6 * 3 * 5 * 23 * 29)^3
    else:
        D = supersingular_D(prime)
        hc_poly = FF['x'](pari(D).polclass())
        root_hc_poly_list = list(hc_poly.roots())

        j_invss = root_hc_poly_list[0][0]
    return FF(j_invss)


@richcmp_method
class SupersingularModule(HeckeModule_free_module):
    r"""
    The module of supersingular points in a given characteristic, with
    given level structure.

    The characteristic must not divide the level.

    .. NOTE:: Currently, only level 1 is implemented.

    EXAMPLES::

        sage: S = SupersingularModule(17)
        sage: S
        Module of supersingular points on X_0(1)/F_17 over Integer Ring
        sage: S = SupersingularModule(16)
        Traceback (most recent call last):
        ...
        ValueError: the argument prime must be a prime number
        sage: S = SupersingularModule(prime=17, level=34)
        Traceback (most recent call last):
        ...
        ValueError: the argument level must be coprime to the argument prime
        sage: S = SupersingularModule(prime=17, level=5)
        Traceback (most recent call last):
        ...
        NotImplementedError: supersingular modules of level > 1 not yet implemented
    """
    def __init__(self, prime=2, level=1, base_ring=ZZ):
        r"""
        Create a supersingular module.

        EXAMPLES::

            sage: SupersingularModule(3)
            Module of supersingular points on X_0(1)/F_3 over Integer Ring
        """
        if not prime.is_prime():
            raise ValueError("the argument prime must be a prime number")
        if prime.divides(level):
            raise ValueError("the argument level must be coprime to the argument prime")
        if level != 1:
            raise NotImplementedError("supersingular modules of level > 1 not yet implemented")
        self.__prime = prime
        self.__finite_field = FiniteField(prime**2, 'a')
        self.__level = level
        self.__hecke_matrices = {}
        HeckeModule_free_module.__init__(self, base_ring,
                                         prime * level, weight=2)

    def _repr_(self) -> str:
        """
        String representation of ``self``.

        EXAMPLES::

            sage: SupersingularModule(11)._repr_()
            'Module of supersingular points on X_0(1)/F_11 over Integer Ring'
        """
        return "Module of supersingular points on X_0(%s)/F_%s over %s" % (
            self.__level, self.__prime, self.base_ring())

    def __richcmp__(self, other, op) -> bool:
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: SupersingularModule(37) == ModularForms(37, 2)
            False
            sage: SupersingularModule(37) == SupersingularModule(37, base_ring=Qp(7))
            False
            sage: SupersingularModule(37) == SupersingularModule(37)
            True
        """
        if not isinstance(other, SupersingularModule):
            return NotImplemented
        return richcmp((self.__level, self.__prime, self.base_ring()),
                       (other.__level, other.__prime, other.base_ring()), op)

    def free_module(self):
        """
        EXAMPLES::

            sage: X = SupersingularModule(37)
            sage: X.free_module()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring

        This illustrates the fix at :issue:`4306`::

            sage: X = SupersingularModule(389)
            sage: X.basis()
            ((1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
             (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
        """
        return ZZ**self.dimension()

    def dimension(self):
        r"""
        Return the dimension of the space of modular forms of weight 2
        and level equal to the level associated to ``self``.

        INPUT:

        - ``self`` -- SupersingularModule object

        OUTPUT: integer; dimension, nonnegative

        EXAMPLES::

            sage: S = SupersingularModule(7)
            sage: S.dimension()
            1

            sage: S = SupersingularModule(15073)
            sage: S.dimension()
            1256

            sage: S = SupersingularModule(83401)
            sage: S.dimension()
            6950

        .. NOTE::

           The case of level > 1 has not yet been implemented.

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        try:
            return self.__dimension
        except AttributeError:
            pass
        if self.__level == 1:
            G = Gamma0(self.__prime)
            self.__dimension = G.dimension_modular_forms(2)
        else:
            raise NotImplementedError
        return self.__dimension

    rank = dimension

    def level(self):
        r"""
        This function returns the level associated to ``self``.

        INPUT:

        - ``self`` -- SupersingularModule object

        OUTPUT: integer; the level, positive

        EXAMPLES::

            sage: S = SupersingularModule(15073)
            sage: S.level()
            1

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        return self.__level

    def prime(self):
        r"""
        Return the characteristic of the finite field associated to ``self``.

        INPUT:

        - ``self`` -- SupersingularModule object

        OUTPUT: integer; characteristic, positive

        EXAMPLES::

            sage: S = SupersingularModule(19)
            sage: S.prime()
            19

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        return self.__prime

    def weight(self):
        r"""
        Return the weight associated to ``self``.

        INPUT:

        - ``self`` -- SupersingularModule object

        OUTPUT: integer; weight, positive

        EXAMPLES::

            sage: S = SupersingularModule(19)
            sage: S.weight()
            2

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        return 2

    def supersingular_points(self):
        r"""
        Compute the supersingular j-invariants over the
        finite field associated to ``self``.

        INPUT:

        - ``self`` -- SupersingularModule object

        OUTPUT:

        - list_j, dict_j -- list_j is the list of supersingular
            j-invariants, dict_j is a dictionary with these
            j-invariants as keys and their indexes as values. The
            latter is used to speed up j-invariant look-up. The
            indexes are based on the order of their *discovery*.

        EXAMPLES:

        The following examples calculate supersingular j-invariants
        over finite fields with characteristic 7, 11 and 37::

            sage: S = SupersingularModule(7)
            sage: S.supersingular_points()
            ([6], {6: 0})

            sage: S = SupersingularModule(11)
            sage: S.supersingular_points()[0]
            [1, 0]

            sage: S = SupersingularModule(37)
            sage: S.supersingular_points()[0]
            [8, 27*a + 23, 10*a + 20]

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        try:
            return (self._ss_points_dic, self._ss_points)
        except AttributeError:
            pass
        Fp2 = self.__finite_field
        level = self.__level
        prime = Fp2.characteristic()
        X = Fp2['x'].gen()
        jinv = supersingular_j(Fp2)

        dim = dimension_supersingular_module(prime, level)

        pos = 0
        # using list to keep track of explored nodes using pos
        ss_points = [jinv]

        # using  to keep track of index of the previous node
        ss_points_pre = [-1]

        # using dictionary for fast j-invariant look-up
        ss_points_dic = {jinv: pos}

        T2_matrix = MatrixSpace(ZZ, dim, sparse=True)(0)

        while pos < len(ss_points):
            if pos == 0:
                neighbors = Phi_polys(2, X, ss_points[pos]).roots()
            else:
                j_prev = ss_points_pre[pos]
                # TODO: These are quadratic polynomials -- maybe we
                # should use the quadratic formula and fast square
                # root finding (??)
                neighbors = Phi2_quad(X, ss_points[j_prev], ss_points[pos]).roots()

            for xj, ej in neighbors:
                if xj not in ss_points_dic:
                    j = len(ss_points)
                    ss_points += [xj]
                    ss_points_pre += [pos]
                    ss_points_dic[xj] = j
                else:
                    j = ss_points_dic[xj]
                T2_matrix[pos, j] += ej
            # end for
            if pos != 0:
                # also record the root from j_prev
                T2_matrix[pos, j_prev] += 1
            pos += 1

        self.__hecke_matrices[2] = T2_matrix
        return (ss_points, ss_points_dic)

    def upper_bound_on_elliptic_factors(self, p=None, ellmax=2):
        r"""
        Return an upper bound (provably correct) on the number of
        elliptic curves of conductor equal to the level of this
        supersingular module.

        INPUT:

        - ``p`` -- (default: 997) prime to work modulo

        ALGORITHM: Currently we only use `T_2`.  Function will be
        extended to use more Hecke operators later.

        The prime p is replaced by the smallest prime that does not
        divide the level.

        EXAMPLES::

            sage: SupersingularModule(37).upper_bound_on_elliptic_factors()
            2

        (There are 4 elliptic curves of conductor 37, but only 2 isogeny
        classes.)
        """
        from sage.misc.verbose import verbose

        # NOTE: The heuristic runtime is *very* roughly `p^2/(2\cdot 10^6)`.
        # ellmax -- (default: 2) use Hecke operators T_ell with ell <= ellmax
        if p is None:
            p = 997

        while self.level() % p == 0:
            p = next_prime(p)

        ell = 2
        t = self.hecke_matrix(ell).change_ring(FiniteField(p))

        # TODO: temporarily try using sparse=False
        # turn this off when sparse rank is optimized.
        t = t.dense_matrix()

        B = ZZ(4 * ell).isqrt()
        bnd = 0
        lower = -B
        upper = B + 1
        for a in range(lower, upper):
            tm = verbose("computing T_%s - %s" % (ell, a))
            if a == lower:
                c = a
            else:
                c = 1
            for i in range(t.nrows()):
                t[i, i] += c
            tm = verbose("computing kernel", tm)
            # dim = t.kernel().dimension()
            dim = t.nullity()
            bnd += dim
            verbose('got dimension = %s; new bound = %s' % (dim, bnd), tm)
        return bnd

    def hecke_matrix(self, L):
        r"""
        Return the `L^{\text{th}}` Hecke matrix.

        INPUT:

        - ``self`` -- SupersingularModule object

        - ``L`` -- integer; positive

        OUTPUT: matrix; sparse integer matrix

        EXAMPLES:

        This example computes the action of the Hecke operator `T_2`
        on the module of supersingular points on `X_0(1)/F_{37}`::

            sage: S = SupersingularModule(37)
            sage: M = S.hecke_matrix(2)
            sage: M
            [1 1 1]
            [1 0 2]
            [1 2 0]

        This example computes the action of the Hecke operator `T_3`
        on the module of supersingular points on `X_0(1)/F_{67}`::

            sage: S = SupersingularModule(67)
            sage: M = S.hecke_matrix(3)
            sage: M
            [0 0 0 0 2 2]
            [0 0 1 1 1 1]
            [0 1 0 2 0 1]
            [0 1 2 0 1 0]
            [1 1 0 1 0 1]
            [1 1 1 0 1 0]

        .. NOTE::

            The first list --- list_j --- returned by the supersingular_points
            function are the rows *and* column indexes of the above hecke
            matrices and its ordering should be kept in mind when interpreting
            these matrices.

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        if L in self.__hecke_matrices:
            return self.__hecke_matrices[L]
        SS, II = self.supersingular_points()
        if L == 2:
            # since T_2 gets computed as a side effect of computing the supersingular points
            return self.__hecke_matrices[2]
        Fp2 = self.__finite_field
        h = len(SS)
        R = self.base_ring()
        T_L = MatrixSpace(R, h)(0)  # mutable
        S, X = Fp2['x'].objgen()

        for i in range(len(SS)):
            ss_i = SS[i]
            phi_L_in_x = Phi_polys(L, X, ss_i)
            rts = phi_L_in_x.roots()
            for r in rts:
                T_L[i, int(II[r[0]])] = r[1]

        self.__hecke_matrices[L] = T_L
        return T_L
