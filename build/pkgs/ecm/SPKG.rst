ecm: Elliptic curve method for integer factorization
====================================================

Description
-----------

GMP-ECM - Elliptic Curve Method for Integer Factorization

Sources can be obtained from https://gitlab.inria.fr/zimmerma/ecm

License
-------

LGPL V3+


Upstream Contact
----------------

-  ecm-discuss@inria.fr

Special Update/Build Instructions
---------------------------------

-  GMP-ECM comes with a self-tuning feature; we could support
   that as an option ($SAGE_TUNE_*=yes) in the future.

-  ECM currently does not (by itself) use the CC and CFLAGS settings
   from 'gmp.h' since we pass (other) options in CFLAGS, and CC is set
   by Sage and might got set by the user. We now at least partially fix
   that
   such that "optimized" code generation options ('-mcpu=...',
   '-mtune=...')
   are used by gcc.
   Of course a user can also manually enable them by setting the
   "global"
   CFLAGS to e.g. '-march=native' on x86[_64] systems, or '-mcpu=...'
   and
   '-mtune=...' on other architectures where "native" isn't supported.
   Note that this doesn't affect the packages' selection of processor-
   specific optimized [assembly] code.
   'spkg-install' already reads the settings from Sage's and also a
   system-wide GMP now, but doesn't (yet) use all of them.
   If SAGE_FAT_BINARY="yes", we should avoid too specific settings of
   "-mcpu=...", and perhaps pass a more generic "--host=..." to
   'configure'.

-  We currently work around a linker bug on MacOS X 10.5 PPC (with
   GCC 4.2.1) which breaks 'configure' if debug symbols are enabled.
   This \*might\* get fixed in later upstream releases.

-  We could save some space by removing the ``src/build.vc10/``
   directory which
   isn't used in Sage. (It gets probably more worth in case also
   directories /
   files for later versions of Microsoft Visual C get added.)
