SAGE_SPKG_CONFIGURE([tox], [
       dnl Early 4.0.x versions have bugs regarding complex factor conditions
       dnl [pkgenv] added in 4.2 - https://tox.wiki/en/latest/upgrading.html#packaging-configuration-and-inheritance
       dnl 4.2.7 for repaired numerical factors
       m4_pushdef([TOX4_MIN_VERSION], [4.2.7])
       AC_CACHE_CHECK([for tox >= ]TOX4_MIN_VERSION, [ac_cv_path_TOX], [
         AC_PATH_PROGS_FEATURE_CHECK([TOX], [tox], [
            tox_version=$($ac_path_TOX --version 2> /dev/null | tail -n 1)
            AX_COMPARE_VERSION([$tox_version], [ge], TOX4_MIN_VERSION, [
                ac_cv_path_TOX="$ac_path_TOX"
                ac_path_TOX_found=:
            ])
         ])
       ])
       AS_IF([test -z "$ac_cv_path_TOX"],
             [sage_spkg_install_tox=yes])
       m4_popdef([TOX4_MIN_VERSION])
])
