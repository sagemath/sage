SAGE_SPKG_CONFIGURE([flint], [
    SAGE_SPKG_DEPCHECK([mpfr], [
        AC_MSG_CHECKING([whether FLINT headers are >= 3.5.0 and < 3.7.0])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
            #include <flint/flint.h>
            #include <flint/padic.h>
            #if __FLINT_RELEASE < 30500 || __FLINT_RELEASE >= 30700
            # error "unsupported FLINT version"
            #endif
        ]])], [
            AC_MSG_RESULT([yes])
            dnl eclib 20261021 uses fmpz_mat_snf_transform, added in FLINT 3.5.
            AC_SEARCH_LIBS([fmpz_mat_snf_transform], [flint], [], [sage_spkg_install_flint=yes])
        ], [
            AC_MSG_RESULT([no])
            sage_spkg_install_flint=yes
        ])
    ])
], [], [], [dnl
     if test x$sage_spkg_install_flint = xyes; then
        AC_SUBST(SAGE_FLINT_PREFIX, ['$SAGE_LOCAL'])
     else
        AC_SUBST(SAGE_FLINT_PREFIX, [''])
     fi
])
