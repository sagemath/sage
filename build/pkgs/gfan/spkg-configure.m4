SAGE_SPKG_CONFIGURE(
    [gfan], [
        AC_CACHE_CHECK([for gfan >= 0.6.2], [ac_cv_path_GFAN], [
        AC_PATH_PROGS_FEATURE_CHECK([GFAN_VERSION], [gfan_version], [
            gfan_version=`$ac_path_GFAN_VERSION | $SED -n '/^gfan/s/gfan//p'`
            AS_IF([test -n "$gfan_version"], [
                AX_COMPARE_VERSION([$gfan_version], [ge], [0.6.2], [
                    ac_cv_path_GFAN_VERSION="$ac_path_GFAN_VERSION"
                    ac_path_GFAN_VERSION_found=:
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_GFAN_VERSION"], [sage_spkg_install_gfan=yes])
])
