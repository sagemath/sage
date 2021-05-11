SAGE_SPKG_CONFIGURE([atlas], [dnl use old test/installation procedure with env. variables
    sage_spkg_install_atlas=yes
    ], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_OPENBLAS])
    AS_IF([test x"$with_blas" = xatlas], [
      sage_require_openblas=no
      sage_require_atlas=yes]) 
])
