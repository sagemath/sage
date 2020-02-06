SAGE_SPKG_CONFIGURE([gsl], [
    m4_pushdef([SAGE_GSL_MINVER],["2.4"])
    SAGE_SPKG_DEPCHECK([atlas openblas], [
      PKG_CHECK_MODULES([GSL], [gsl >= $SAGE_GSL_MINVER], [
        PKG_CHECK_VAR([GSLPCDIR], [gsl], [pcfiledir], [
          AC_CONFIG_COMMANDS([GSLPCPROCESS], [
            $SED -e 's/\${GSL_CBLAS_LIB}\ //' \
                 -e 's/GSL_CBLAS_LIB.*/Requires: cblas/' $GSL_PC \
                  > "$SAGE_SRC"/lib/pkgconfig/gsl.pc
          ], [
            SED=$ac_cv_path_SED
            GSL_PC="$GSLPCDIR"/gsl.pc
          ])
        ], [
        AC_MSG_WARN([Unable to locate the directory of gsl.pc. This should not happen!])
       sage_spkg_install_gsl=yes
       ])
      ], [sage_spkg_install_gsl=yes])
    ])
    m4_popdef([SAGE_GSL_MINVER])
])
