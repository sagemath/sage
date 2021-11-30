SAGE_SPKG_CONFIGURE([gsl], [
    m4_pushdef([SAGE_GSL_MINVER],["2.4"])
    SAGE_SPKG_DEPCHECK([openblas], [
      PKG_CHECK_MODULES([GSL], [gsl >= $SAGE_GSL_MINVER], [
        PKG_CHECK_VAR([GSLPCDIR], [gsl], [pcfiledir], [
          GSL_PC="$GSLPCDIR"/gsl.pc
          AC_SUBST([SAGE_SYSTEM_FACADE_PC_FILES])
          AS_VAR_APPEND([SAGE_SYSTEM_FACADE_PC_FILES], [" \$(SAGE_PKGCONFIG)/gsl.pc"])
          AC_SUBST([SAGE_GSL_PC_COMMAND],["\$(SED) -e 's/\$\${GSL_CBLAS_LIB}//' -e \"s/^GSL_CBLAS_LIB=.*/Requires: cblas/\" \"$GSL_PC\" > \"\$(@)\""])
        ], [
        AC_MSG_WARN([Unable to locate the directory of gsl.pc. This should not happen!])
       sage_spkg_install_gsl=yes
       ])
      ], [sage_spkg_install_gsl=yes])
    ])
    m4_popdef([SAGE_GSL_MINVER])
])
