SAGE_SPKG_CONFIGURE([eclib], [
  SAGE_SPKG_DEPCHECK([ntl pari flint], [
    dnl use existing eclib only if the version reported by pkg-config is recent enough
    m4_pushdef([SAGE_ECLIB_VER],["20231212"])
    PKG_CHECK_MODULES([ECLIB], [eclib >= SAGE_ECLIB_VER], [
      AC_CACHE_CHECK([for mwrank version == SAGE_ECLIB_VER], [ac_cv_path_MWRANK], [
        AC_PATH_PROGS_FEATURE_CHECK([MWRANK], [mwrank], [
            mwrank_version=`$ac_path_MWRANK -V 2>&1`
            AX_COMPARE_VERSION([$mwrank_version], [ge], [SAGE_ECLIB_VER], [
                ac_cv_path_MWRANK="$ac_path_MWRANK"
            ])
        ])
      ])
      AS_IF([test -z "$ac_cv_path_MWRANK"], [sage_spkg_install_eclib=yes])
    ], [
    sage_spkg_install_eclib=yes])
  ])
  m4_popdef([SAGE_ECLIB_VER])
])
