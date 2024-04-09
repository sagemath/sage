Fixed and Arbitrary Precision Numerical Fields
==============================================

Floating-Point Arithmetic
-------------------------

Sage supports arbitrary precision real (:class:`RealField`) and complex fields
(:class:`ComplexField`). Sage also provides two optimized fixed precision fields for
numerical computation, the real double (:class:`RealDoubleField`) and complex double
fields (:class:`ComplexDoubleField`).

Real and complex double elements are optimized implementations that use the
:ref:`GNU Scientific Library <spkg_gsl>` for arithmetic and some special functions.
Arbitrary precision real and complex numbers are implemented using the
:ref:`MPFR <spkg_mpfr>` library, which builds on :ref:`GMP <spkg_gmp>`.
In many cases, the :ref:`PARI <spkg_pari>` C-library is used to compute
special functions when implementations aren't otherwise available.

.. toctree::
   :maxdepth: 1

   sage/rings/real_mpfr
   sage/rings/complex_mpfr
   sage/rings/complex_mpc
   sage/rings/real_double
   sage/rings/complex_double

Interval Arithmetic
-------------------

Sage implements real and complex interval arithmetic using
:ref:`MPFI <spkg_mpfi>` (:class:`RealIntervalField`, :class:`ComplexIntervalField`)
and :ref:`FLINT <spkg_flint>` (:class:`RealBallField`, :class:`ComplexBallField`).

.. toctree::
   :maxdepth: 1

   sage/rings/real_mpfi
   sage/rings/real_interval_absolute
   sage/rings/complex_interval_field
   sage/rings/complex_interval

   sage/rings/real_arb
   sage/rings/complex_arb

Exact Real Arithmetic
---------------------

.. toctree::
   :maxdepth: 1

   sage/rings/real_lazy

.. include:: ../footer.txt

