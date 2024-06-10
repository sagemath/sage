mpfr: Multiple-precision floating-point computations with correct rounding
==========================================================================

Description
-----------

The MPFR library is a C library for multiple-precision floating-point
computations with correct rounding. MPFR has continuously been supported
by the INRIA and the current main authors come from the Caramba and AriC
project-teams at Loria (Nancy, France) and LIP (Lyon, France)
respectively; see more on the credit page. MPFR is based on the GMP
multiple-precision library.

The main goal of MPFR is to provide a library for multiple-precision
floating-point computation which is both efficient and has a
well-defined semantics. It copies the good ideas from the ANSI/IEEE-754
standard for double-precision floating-point arithmetic (53-bit
significand).

License
-------

MPFR is free. It is distributed under the GNU Lesser General Public
License (GNU Lesser GPL), version 3 or later (2.1 or later for MPFR
versions until 2.4.x). The library has been registered in France by the
Agence de Protection des Programmes under the number IDDN FR 001 120020
00 R P 2000 000 10800, on 15 March 2000. This license guarantees your
freedom to share and change MPFR, to make sure MPFR is free for all its
users. Unlike the ordinary General Public License, the Lesser GPL
enables developers of non-free programs to use MPFR in their programs.
If you have written a new function for MPFR or improved an existing one,
please share your work!


Upstream Contact
----------------

The MPFR website is located at http://mpfr.org/

The MPFR team can be contacted via the MPFR mailing list: mpfr@loria.fr

Special Update/Build Instructions
---------------------------------

-  Make sure MPFR's settings of ``CC`` and ``CFLAGS`` still get properly
   extracted,
   currently from its ``config.log`` in the ``src/`` directory.

-  We should remove the ``configure`` option ``--disable-thread-safe``
   in case
   the issues without that have meanwhile been fixed. (Then we should
   actually pass ``--enable-thread-safe``.)

TODO
----

-  ``--disable-thread-safe`` should be switched to ``--enable-thread-safe``,
   need to check that this works on the buildbot machines
