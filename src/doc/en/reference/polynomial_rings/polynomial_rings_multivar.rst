
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

are implemented using the PolyBoRi library (cf. :mod:`sage.rings.polynomial.pbori`).


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
