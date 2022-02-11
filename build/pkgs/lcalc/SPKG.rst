lcalc: L-function calculator
============================

Description
-----------

Michael Rubinstein's L-function calculator.

License
-------

-  LGPL V2+


Upstream contact
----------------

Michael Rubinstein <mrubinst@uwaterloo.ca>

Sources: http://oto.math.uwaterloo.ca/~mrubinst/L_function_public/L.html

Newer beta version 1.3 (not yet in Sage):
http://code.google.com/p/l-calc/

Dependencies
------------

-  GMP/MPIR
-  MPFR
-  PARI
-  GNU patch


Special Update/Build Instructions
---------------------------------

-  There is some garbage in the upstream sources which should be
   removed::

     src/include/.Lexplicit_formula.h.swp
     src/include/.Lvalue.h.swp
     src/include/._.DS_Store
     src/include/.DS_Store
     src/include/Lexplicit_formula.h.swap.crap
     src/include/Lvalue.h.bak
     src/src/Makefile.old
     src/src/.Makefile.old.swp
     src/src/._.DS_Store
     src/src/.DS_Store
     src/src/.Lcommandline.ggo.swp
     src/src/libLfunction.a

-  We (and apparently also upstream) currently don't build Lcalc's tests
   (see Makefile), hence there's no spkg-check.
   This might change in newer upstream versions.

-  The original Makefile uses $(CC) to compile C++ (also using
   $(CCFLAGS)),
   which it defines to 'g++', and hardcodes 'g++' when linking the
   shared
   library. (It should use $(CXX) instead, which might \*default\* to
   'g++'.)
   We now (lcalc-1.23.p10) patch the Makefile also to use $(CXX) for
   compiling
   and linking C++; $(CXX) now \*defaults\* to 'g++', and $(CC) to
   'gcc', but
   both can be overridden by simply setting their respective environment
   variables. (Same for $(INSTALL_DIR) btw.)

Patches
-------

-  Makefile.patch:

   We change a lot there, since Lcalc doesn't have a 'configure' script,
   and hence the Makefile is supposed to be edited to customize Lcalc
   (build
   options, locations of headers and libraries etc.).
   Besides that, we

   -  put CXXFLAGS into Lcalc's "CCFLAGS" used for compiling C++,
   -  remove some stuff involving LDFLAGS1 and LDFLAGS2, setting just
      LDFLAGS,
   -  use $(MAKE) instead of 'make' in the crude build receipts,
   -  use CXXFLAG64 when linking the shared library,
   -  now use $(CXX) for compiling and linking C++, which \*defaults\* to
      'g++',
      but can be overridden by setting the environment variable of the same
      name. ($(CC) now \*defaults\* to 'gcc', although currently not really
      used as far as I can see.)
   -  $(INSTALL_DIR) can now be overridden by simply setting the
      environment
      variable of the same name.

-  Lcommon.h.patch:

   Uncomment the definition of lcalc_to_double(const long double& x).
   (Necessary for GCC >= 4.6.0, cf. #10892.)
   Comment from there:
   The reason is the following code horror from
   src/src/include/Lcommon.h:
   [...]
   But somebody who is familiar with the codebase should really rewrite
   lcalc
   to not redefine the double() cast, thats just fragile and will sooner
   or
   later again fail inside some system headers.

-  pari-2.7.patch:

   Various changes to port to newer versions of PARI.

-  time.h.patch:

   (Patches src/include/Lcommandline_numbertheory.h)
   Include also <time.h> in Lcommandline_numbertheory.h (at least
   required
   on Cygwin, cf. #9845).
   This should get reported upstream.

-  lcalc-1.23_default_parameters_1.patch: Make Lcalc (1.23) build with
   GCC 4.9
