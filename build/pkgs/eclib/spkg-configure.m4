SAGE_SPKG_CONFIGURE([eclib], [
  SAGE_SPKG_DEPCHECK([ntl pari flint], [
    PKG_CHECK_EXISTS([eclib = 20261021], [eclib_version=20261021], [eclib_version=20250627])
    PKG_CHECK_MODULES([ECLIB], [eclib = $eclib_version], [
      AC_CACHE_CHECK([for mwrank version >= 20250627], [ac_cv_path_MWRANK], [
        AC_PATH_PROGS_FEATURE_CHECK([MWRANK], [mwrank], [
            mwrank_version=`$ac_path_MWRANK -V 2>&1`
            AX_COMPARE_VERSION([$mwrank_version], [ge], [20250627], [
                ac_cv_path_MWRANK="$ac_path_MWRANK"
            ])
        ])
      ])
      AS_IF([test -z "$ac_cv_path_MWRANK"], [sage_spkg_install_eclib=yes])
    ], [
      sage_spkg_install_eclib=yes])
  ])
])
