SAGE_SPKG_CONFIGURE(
    [ninja_build], [
        AC_CACHE_CHECK([for ninja >= 1.7.2], [ac_cv_path_NINJA], [
        AC_PATH_PROGS_FEATURE_CHECK([NINJA], [ninja], [
            ninja_version=`$ac_path_NINJA --version 2>&1 \
                | $SED -n -e 's/\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\).*/\1/p'`
            AS_IF([test -n "$ninja_version"], [
                AX_COMPARE_VERSION([$ninja_version], [ge], [1.7.2], [
                    ac_cv_path_NINJA="$ac_path_NINJA"
                    ac_path_NINJA_found=:
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_NINJA"], [sage_spkg_install_ninja_build=yes])
])
