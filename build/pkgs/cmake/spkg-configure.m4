SAGE_SPKG_CONFIGURE(
    [cmake], [
	AC_CACHE_CHECK([for cmake >= 3.4], [ac_cv_path_CMAKE], [
        AC_PATH_PROGS_FEATURE_CHECK([CMAKE], [cmake], [
            cmake_version=`$ac_path_CMAKE --version 2>&1 \
                | $SED -n -e 's/cmake version *\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\)/\1/p'`
            AS_IF([test -n "$cmake_version"], [
                AX_COMPARE_VERSION([$cmake_version], [ge], [3.4], [
                    ac_cv_path_CMAKE="$ac_path_CMAKE"
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_CMAKE"], [sage_spkg_install_cmake=yes])
])
