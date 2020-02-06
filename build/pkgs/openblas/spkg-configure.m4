SAGE_SPKG_CONFIGURE([openblas], [
 SAGE_SPKG_DEPCHECK([gfortran], [
  PKG_CHECK_MODULES([OPENBLAS], [openblas >= 0.2.20], [
    PKG_CHECK_VAR([OPENBLASPCDIR], [openblas], [pcfiledir], [
       sage_install_blas_pc=yes
       AC_SEARCH_LIBS([cblas_dgemm], [openblas], [dnl openblas works as cblas
             sage_install_cblas_pc=yes
             ], [
             dnl openblas does not work as cblas; try to use system's cblas as is
             PKG_CHECK_MODULES([CBLAS], [cblas], [], [sage_spkg_install_openblas=yes])
          ])
       AC_FC_FREEFORM([AC_FC_FUNC([dgeqrf])])
       AC_SEARCH_LIBS([$dgeqrf], [openblas], [dnl openblas works as lapack
             sage_install_lapack_pc=yes
             ], [
             dnl openblas does not work as lapack; try to use system's lapack as is
             PKG_CHECK_MODULES([LAPACK], [lapack], [], [sage_spkg_install_openblas=yes])
          ])
       ], [
       AC_MSG_WARN([Unable to locate the directory of openblas.pc. This should not happen!])
       sage_spkg_install_openblas=yes
       ])
    ], [sage_spkg_install_openblas=yes])
    AS_IF([test x$sage_spkg_install_openblas != xyes], [
       m4_foreach([blaslibnam], [blas, cblas, lapack], [
        AS_IF([test x$sage_install_]blaslibnam[_pc = xyes], [
         AC_CONFIG_LINKS([$SAGE_SRC/lib/pkgconfig/]blaslibnam[.pc:$OPENBLASPCDIR/openblas.pc])])
       ])
    ])
 ])
 ], [
  AS_IF([test "x$with_blas" = xopenblas], [
     sage_require_openblas=yes
     sage_require_atlas=no])
  ], [
  AC_MSG_CHECKING([BLAS library])
  AC_ARG_WITH([blas],
  [AS_HELP_STRING([--with-blas=openblas],
    [use OpenBLAS as BLAS library (default)])]
  [AS_HELP_STRING([--with-blas=atlas],
    [use ATLAS as BLAS library])],,
    [with_blas=openblas]  # default
  )
  AS_CASE(["$with_blas"],
    [openblas], [],
    [atlas],    [sage_spkg_install_openblas=no],
                [AC_MSG_ERROR([allowed values for --with-blas are 'atlas' and 'openblas'])])
  AC_MSG_RESULT([$with_blas])
  AC_SUBST([SAGE_BLAS], [$with_blas])
  ]
)
