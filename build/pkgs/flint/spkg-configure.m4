SAGE_SPKG_CONFIGURE([flint], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_MPFR])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_NTL])
    AC_MSG_CHECKING([installing mpfr or ntl? ])
    if test x$sage_spkg_install_mpfr = xyes -o x$sage_spkg_install_ntl = xyes; then
        AC_MSG_RESULT([yes; install flint as well])
        sage_spkg_install_flint=yes
    else
        AC_CHECK_HEADER(flint/flint.h, [
          dnl fmpz_mat_is_hadamard appears in Flint 2.5.0
          AC_SEARCH_LIBS([fmpz_mat_is_hadamard], [flint], [
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
		         sage_spkg_install_flint=yes])
            ], [sage_spkg_install_flint=yes])
          ], [sage_spkg_install_flint=yes])
        ], [sage_spkg_install_flint=yes])
    fi
], [], [], [
     if test x$sage_spkg_install_flint = xyes; then
        AC_SUBST(SAGE_FLINT_PREFIX, ['$SAGE_LOCAL'])
        AC_MSG_RESULT([using Sage's flint SPKG])
     else
        AC_SUBST(SAGE_FLINT_PREFIX, [''])
        AC_MSG_RESULT([using flint library from the system])
     fi
])
