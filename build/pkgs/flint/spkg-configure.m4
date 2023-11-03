SAGE_SPKG_CONFIGURE([flint], [
    SAGE_SPKG_DEPCHECK([mpfr ntl], [
        AC_CHECK_HEADER(flint/flint.h, [
          dnl flint_parallel_binary_splitting appears in Flint 2.9.0, needed by arb 2.23
          AC_SEARCH_LIBS([flint_parallel_binary_splitting], [flint], [
            dnl check that NTL is linked in
            AC_SEARCH_LIBS([fmpz_poly_get_ZZX], [flint], [

              AC_MSG_CHECKING([that GC is not enabled in Flint... ])
              AC_RUN_IFELSE([
                 AC_LANG_PROGRAM([[#include <flint/flint.h>]], [
                                  [#ifdef HAVE_GC]
                                     [return HAVE_GC;]
                                  [#else]
                                     [return 0;]
                                  [#endif]])],
                 [AC_MSG_RESULT([GC not enabled. Good.])],
                        [AC_MSG_RESULT([GC enabled. Incompatible with Sage.])
                         sage_spkg_install_flint=yes],
                 [AC_MSG_RESULT(["cross compiling. assuming GC is not enabled"])])
            ], [sage_spkg_install_flint=yes])
          ], [sage_spkg_install_flint=yes])
        ], [sage_spkg_install_flint=yes])
    ])
], [], [], [
     if test x$sage_spkg_install_flint = xyes; then
        AC_SUBST(SAGE_FLINT_PREFIX, ['$SAGE_LOCAL'])
     else
        AC_SUBST(SAGE_FLINT_PREFIX, [''])
     fi
])
