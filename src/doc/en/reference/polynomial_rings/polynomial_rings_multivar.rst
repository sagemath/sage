
Multivariate Polynomials and Polynomial Rings
=============================================

Sage implements multivariate polynomial rings through several
backends. The most generic implementation uses the classes
:class:`sage.rings.polynomial.polydict.PolyDict` and
:class:`sage.rings.polynomial.polydict.ETuple` to construct a dictionary with
exponent tuples as keys and coefficients as values.

Additionally, specialized and optimized implementations over many
specific coefficient rings are implemented via a shared library interface to
SINGULAR; and polynomials in the boolean polynomial ring

.. math::

    \GF{2}[x_1,...,x_n]/ \langle x_1^2+x_1,...,x_n^2+x_n \rangle.

are implemented using the PolyBoRi library (cf.
:mod:`sage.rings.polynomial.pbori.pbori`).


.. toctree::
   :maxdepth: 1

   sage/rings/polynomial/term_order

   sage/rings/polynomial/multi_polynomial_ring_base
   sage/rings/polynomial/multi_polynomial

   sage/rings/polynomial/multi_polynomial_ring
   sage/rings/polynomial/multi_polynomial_element
   sage/rings/polynomial/multi_polynomial_ideal

   sage/rings/polynomial/multi_polynomial_sequence

   sage/rings/polynomial/multi_polynomial_libsingular
   sage/rings/polynomial/multi_polynomial_ideal_libsingular

   sage/rings/polynomial/msolve

   sage/rings/polynomial/polydict
   sage/rings/polynomial/hilbert

   sage/rings/polynomial/flatten

   sage/rings/monomials


Absolute factorization and absolute irreducibility
==================================================

In computational algebra and algebraic geometry, **absolute factorization**
refers to factorization of polynomials over the algebraic closure
:math:`\overline{\QQ}` of the rational numbers. A polynomial is called
**absolutely irreducible** if it is irreducible over
:math:`\overline{\QQ}`.

SageMath supports absolute factorization natively via the algebraic number
fields ``QQbar`` (the field of algebraic numbers) and ``AA`` (the field of
real algebraic numbers).

Multivariate absolute factorization over ``QQbar``
--------------------------------------------------

For multivariate polynomials over ``QQbar``, Sage computes an absolute
factorization using Singular’s ``absfact`` library::

    sage: R.<x,y> = QQbar[]
    sage: f = x^2 + y^2
    sage: f.factor()
    (x - I*y) * (x + I*y)

This shows that ``x^2 + y^2`` is reducible over
:math:`\overline{\QQ}`, even though it is irreducible over ``QQ``.

Absolute irreducibility test
----------------------------

A common way to test whether a polynomial is absolutely irreducible is
to factor it after extending the coefficient ring to ``QQbar``::

    sage: R.<x,y> = QQ[]
    sage: f = x^2 + y^2
    sage: f.is_irreducible()
    True
    sage: f.change_ring(QQbar).is_irreducible()
    False

Thus, ``f`` is irreducible over ``QQ`` but **not absolutely irreducible**.

Real algebraic coefficients
---------------------------

When factoring over ``AA`` (the field of real algebraic numbers), complex
conjugate factors are combined whenever possible, so that the result
lies in ``AA``::

    sage: R.<x,y> = AA[]
    sage: (x^2 + y^2).factor()
    x^2 + y^2
