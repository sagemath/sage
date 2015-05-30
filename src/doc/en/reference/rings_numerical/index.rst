Fixed and Arbitrary Precision Numerical Fields
==============================================

Sage supports two optimized fixed precision fields for numerical
computation, the real double (RealDoubleField) and complex double
fields (ComplexDoubleField).
Sage also supports arbitrary precision real (RealField) and
complex fields (ComplexField), and real and complex
interval arithmetic (RealIntervalField and ComplexIntervalField).


Real and complex double elements are optimized implementations that
use the GNU Scientific Library for arithmetic and some special
functions.  Arbitrary precision real and complex numbers are
implemented using the MPFR library, which builds on GMP.  (Note that
Sage doesn't currently use the MPC library.)  The interval arithmetic
field is implemented using the MPFI library.

In many cases the PARI C-library is used to compute special functions
when implementations aren't otherwise available.

.. toctree::
   :maxdepth: 2

   sage/rings/real_double
   sage/rings/real_mpfr
   sage/rings/real_mpfi
   sage/rings/real_interval_field
   sage/rings/real_interval_absolute
   sage/rings/real_lazy

   sage/rings/complex_double
   sage/rings/complex_field
   sage/rings/complex_number
   sage/rings/complex_mpc
   sage/rings/complex_interval_field
   sage/rings/complex_interval

.. Modules depending on optional packages:
.. sage/rings/real_arb

.. include:: ../footer.txt

