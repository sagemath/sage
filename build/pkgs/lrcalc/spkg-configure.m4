SAGE_SPKG_CONFIGURE([lrcalc], [
    AC_LANG_PUSH([C])
    AC_CHECK_HEADERS([lrcalc/schublib.h] [lrcalc/lrcoef.h], [
          AC_SEARCH_LIBS([mult_poly_schubert], [lrcalc], [], [sage_spkg_install_lrcalc=yes])
    ], [sage_spkg_install_lrcalc=yes],
    [[#ifndef _IVECTOR_H
      #include <lrcalc/ivector.h>
      #endif]])
    AC_LANG_POP
])
