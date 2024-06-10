SAGE_SPKG_CONFIGURE([iml], [
  SAGE_SPKG_DEPCHECK([gmp openblas], [
    AC_CHECK_HEADER([iml.h], [
      AC_SEARCH_LIBS([nonsingSolvLlhsMM], [iml], [],
                     [sage_spkg_install_iml=yes])
    ], [
      sage_spkg_install_iml=yes
    ], [
      #ifdef HAVE_GMP_H
      #include <gmp.h>
      #endif
    ])
  ])
])
