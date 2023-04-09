SAGE_SPKG_CONFIGURE(
    [meson], [
        AC_CACHE_CHECK([for meson >= 0.64], [ac_cv_path_MESON], [
        AC_PATH_PROGS_FEATURE_CHECK([MESON], [meson], [
            meson_version=`$ac_path_MESON --version 2>&1`
            AS_IF([test -n "$meson_version"], [
                AX_COMPARE_VERSION([$meson_version], [ge], [0.64], [
                    ac_cv_path_MESON="$ac_path_MESON"
                    ac_path_MESON_found=:
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_MESON"], [sage_spkg_install_meson=yes])
])
