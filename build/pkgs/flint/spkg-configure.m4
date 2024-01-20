SAGE_SPKG_CONFIGURE([flint], [
    SAGE_SPKG_DEPCHECK([mpfr], [
        AC_CHECK_HEADER(flint/flint.h, [
          dnl gr_get_fexpr appears in Flint 3.0
          AC_SEARCH_LIBS([gr_get_fexpr], [flint], [], [sage_spkg_install_flint=yes])
        ], [sage_spkg_install_flint=yes])
    ])
], [], [], [
     if test x$sage_spkg_install_flint = xyes; then
        AC_SUBST(SAGE_FLINT_PREFIX, ['$SAGE_LOCAL'])
     else
        AC_SUBST(SAGE_FLINT_PREFIX, [''])
     fi
])
