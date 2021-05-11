Fixed and Arbitrary Precision Numerical Fields
==============================================

Floating-Point Arithmetic
-------------------------

Sage supports arbitrary precision real (RealField) and complex fields
(ComplexField). Sage also provides two optimized fixed precision fields for
numerical computation, the real double (RealDoubleField) and complex double
fields (ComplexDoubleField).

Real and complex double elements are optimized implementations that use the
GNU Scientific Library for arithmetic and some special functions.  Arbitrary
precision real and complex numbers are implemented using the MPFR library,
which builds on GMP. In many cases the PARI C-library is used to compute
special functions when implementations aren't otherwise available.

.. toctree::
   :maxdepth: 2

   sage/rings/real_mpfr
   sage/rings/complex_mpfr
   sage/rings/complex_mpc
   sage/rings/real_double
   sage/rings/complex_double

Interval Arithmetic
-------------------

Sage implements real and complex interval arithmetic using MPFI
(RealIntervalField, ComplexIntervalField) and arb (RealBallField,
ComplexBallField).

.. toctree::
   :maxdepth: 2

   sage/rings/real_mpfi
   sage/rings/real_interval_absolute
   sage/rings/complex_interval_field
   sage/rings/complex_interval

   sage/rings/real_arb
   sage/rings/complex_arb

Exact Real Arithmetic
---------------------

.. toctree::
   :maxdepth: 2

   sage/rings/real_lazy

.. include:: ../footer.txt

