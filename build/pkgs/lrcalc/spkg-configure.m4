SAGE_SPKG_CONFIGURE([lrcalc], [
    AC_CHECK_HEADERS([lrcalc/schublib.h], [
          AC_SEARCH_LIBS([mult_poly_schubert], [lrcalc], [], [sage_spkg_install_lrcalc=yes])
    ], [sage_spkg_install_lrcalc=yes],
    [[#ifndef _HASHTAB_H
      #include <lrcalc/hashtab.h>
      #include <lrcalc/vector.h>
      #endif]])
])
