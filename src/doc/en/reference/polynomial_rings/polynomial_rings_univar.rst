
Univariate Polynomials and Polynomial Rings
===========================================

Sage's architecture for polynomials 'under the hood' is complex, interfacing to
a variety of C/C++ libraries for polynomials over specific rings. In practice,
the user rarely has to worry about which backend is being used.

The hierarchy of class inheritance is somewhat confusing, since most of the
polynomial element classes are implemented as Cython extension types rather
than pure Python classes and thus can only inherit from a single base class,
whereas others have multiple bases.

.. toctree::
   :maxdepth: 2

   sage/rings/polynomial/polynomial_ring
   sage/rings/polynomial/polynomial_ring_homomorphism

   sage/rings/polynomial/polynomial_element
   sage/rings/polynomial/polynomial_element_generic
   sage/rings/polynomial/polynomial_gf2x
   sage/rings/polynomial/polynomial_number_field
   sage/rings/polynomial/polynomial_integer_dense_flint
   sage/rings/polynomial/polynomial_integer_dense_ntl
   sage/rings/polynomial/polynomial_rational_flint
   sage/rings/polynomial/polynomial_zmod_flint
   sage/rings/polynomial/polynomial_modn_dense_ntl
   sage/rings/polynomial/polynomial_real_mpfr_dense
   sage/rings/polynomial/polynomial_singular_interface
   sage/rings/polynomial/padics/polynomial_padic
   sage/rings/polynomial/padics/polynomial_padic_capped_relative_dense
   sage/rings/polynomial/padics/polynomial_padic_flat
   sage/rings/polynomial/polynomial_zz_pex

   sage/rings/polynomial/real_roots
   sage/rings/polynomial/complex_roots
   sage/rings/polynomial/refine_root

   sage/rings/polynomial/ideal
   sage/rings/polynomial/polynomial_quotient_ring
   sage/rings/polynomial/polynomial_quotient_ring_element

   sage/rings/polynomial/polynomial_compiled
   sage/rings/polynomial/polynomial_fateman
