mpfi: Multiple precision interval arithmetic library based on MPFR
==================================================================

Description
-----------

MPFI is a library for interval arithmetic, which is built upon the MPFR
multiple precision floating-point arithmetic.

MPFI is intended to be a portable library written in C for arbitrary
precision interval arithmetic with intervals represented using MPFR
reliable floating-point numbers. It is based on the GNU MP library and
on the MPFR library. The purpose of an arbitrary precision interval
arithmetic is on the one hand to get "guaranteed" results, thanks to
interval computation, and on the other hand to obtain accurate results,
thanks to multiple precision arithmetic. The MPFI library is built upon
MPFR in order to benefit from the correct rounding provided, for each
operation or function, by MPFR. Further advantages of using MPFR are its
portability and compliance with the IEEE 754 standard for floating-point
arithmetic.

License
-------

This version of MPFI is released under the GNU Lesser General Public
License. It is permitted to link MPFI to non-free programs, as long as
when distributing them the MPFI source code and a means to re-link with
a modified MPFI is provided.


Upstream Contact
----------------

http://perso.ens-lyon.fr/nathalie.revol/software.html

The MPFI website is located at https://gitlab.inria.fr/mpfi/mpfi

The MPFI team can be contacted via the MPFI mailing list: mpfi-users@inria.fr

