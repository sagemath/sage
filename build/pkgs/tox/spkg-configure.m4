SAGE_SPKG_CONFIGURE([tox], [
       dnl Use non-ancient tox with full support for PEP 517.
       m4_pushdef([TOX_MIN_VERSION], [3.21.4])
       AC_CACHE_CHECK([for tox >= ]TOX_MIN_VERSION, [ac_cv_path_TOX], [
         AC_PATH_PROGS_FEATURE_CHECK([TOX], [tox], [
            tox_version=$($ac_path_TOX --version 2> /dev/null | tail -1)
            AS_IF([test -n "$tox_version"], [
                AX_COMPARE_VERSION([$tox_version], [ge], TOX_MIN_VERSION, [
                    ac_cv_path_TOX="$ac_path_TOX"
                    ac_path_TOX_found=:
                ])
            ])
         ])
       ])
       AS_IF([test -z "$ac_cv_path_TOX"],
             [sage_spkg_install_tox=yes])
])
