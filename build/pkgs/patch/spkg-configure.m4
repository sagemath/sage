SAGE_SPKG_CONFIGURE(
    [patch], [
        AC_CACHE_CHECK([for GNU patch >= 2.7.0], [ac_cv_path_PATCH], [
        AC_PATH_PROGS_FEATURE_CHECK([PATCH], [patch], [
            patch_version=`$ac_path_PATCH --version 2>&1 \
                | $SED -n -e 's/GNU patch *\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\)/\1/p'`
            AS_IF([test -n "$patch_version"], [
                AX_COMPARE_VERSION([$patch_version], [ge], [2.7.0], [
                    ac_cv_path_PATCH="$ac_path_PATCH"
                    ac_path_PATCH_found=:
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_PATCH"], [sage_spkg_install_patch=yes])
])
