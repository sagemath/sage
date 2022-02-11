from .PyPolyBoRi import Polynomial, BoolePolynomialVector
from .pbori import FGLMStrategy, BooleSet


def _fglm(I, from_ring, to_ring):
    r"""
    Unchecked variant of fglm
    """
    vec = BoolePolynomialVector(I)
    return FGLMStrategy(from_ring, to_ring, vec).main()


def fglm(I, from_ring, to_ring):
    r"""
    Convert *reduced* Groebner Basis in ``from_ring`` to a GroebnerBasis
    in ``to_ring``.

    It acts independent of the global ring, which is restored at the end of the
    computation.

    TESTS::

        sage: from sage.rings.polynomial.pbori import *
        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import OrderCode
        sage: dp_asc = OrderCode.dp_asc
        sage: r=declare_ring(['x','y','z'],dict())
        sage: old_ring = r
        sage: new_ring = old_ring.clone(ordering=dp_asc)
        sage: (x,y,z) = [old_ring.variable(i) for i in range(3)]
        sage: ideal=[x+z, y+z]# lp Groebner basis
        sage: from sage.rings.polynomial.pbori.fglm import fglm
        sage: list(fglm(ideal, old_ring, new_ring))
        [y + x, z + x]
    """
    for poly in I:
        if poly.ring().id() != from_ring.id():
            raise ValueError("Ideal I must be from the first ring argument")
    return _fglm(I, from_ring, to_ring)


def vars_real_divisors(monomial, monomial_set):
    r"""
    Return all elements of ``monomial_set``, which result multiplied by a variable in monomial.

    EXAMPLES::

        sage: from sage.rings.polynomial.pbori.pbori import *
        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import OrderCode
        sage: dp_asc = OrderCode.dp_asc
        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import Ring
        sage: r = Ring(1000)
        sage: x = r.variable
        sage: b = BooleSet([x(1)*x(2),x(2)])
        sage: from sage.rings.polynomial.pbori.fglm import vars_real_divisors
        sage: vars_real_divisors(x(1)*x(2)*x(3),b)
        {{x(1),x(2)}}
    """
    return BooleSet(Polynomial(monomial_set.divisors_of(monomial)).
                    graded_part(monomial.deg() - 1))


def m_k_plus_one(completed_elements, variables):
    r"""
    Calculate $m_{k+1}$ from the FGLM algorithm.

    Calculate $m_{k+1}$ from the FGLM algorithm as described in Wichmann [Wich1997]_.

    .. NOTE::

        It would be nice to be able to efficiently extract the smallest term of a polynomial.

    EXAMPLES::

        sage: from sage.rings.polynomial.pbori.pbori import *
        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import OrderCode
        sage: from sage.rings.polynomial.pbori.fglm import m_k_plus_one
        sage: dp_asc = OrderCode.dp_asc
        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import Ring
        sage: r = Ring(1000)
        sage: x = r.variable
        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import Monomial
        sage: s = BooleSet([x(1)*x(2),x(1),x(2),Monomial(r),x(3)])
        sage: variables = BooleSet([x(1),x(2),x(3)])
        sage: m_k_plus_one(s,variables)
        x(2)*x(3)
        sage: r2 = r.clone(ordering=dp_asc)
        sage: m_k_plus_one(r2(s).set(),r2(variables).set())
        x(1)*x(3)
    """
    return sorted(completed_elements.cartesian_product(variables).diff(
        completed_elements))[0]
