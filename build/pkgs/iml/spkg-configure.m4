SAGE_SPKG_CONFIGURE([iml], [
    dnl m4_pushdef([SAGE_IML_MINVER],["1.0.4"]) 1.0.5 seems OK too
    SAGE_SPKG_DEPCHECK([gmp mpir openblas], [
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
    dnl m4_popdef([SAGE_IML_MINVER])
])
