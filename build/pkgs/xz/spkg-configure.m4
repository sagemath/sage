SAGE_SPKG_CONFIGURE([xz], [
    AC_CHECK_LIB([lzma], [lzma_raw_decoder], [lzma_cv_liblzma=yes], [lzma_cv_liblzma=no])
    AC_CHECK_HEADER([lzma.h], [lzma_cv_lzma_h=yes], [lzma_cv_lzma_h=no])
    if test "$lzma_cv_liblzma" = "yes" && test "$lzma_cv_lzma_h" = "yes"; then
        AC_CACHE_CHECK([for xz >= 4.999.0], [ac_cv_path_XZ], [
            AC_PATH_PROGS_FEATURE_CHECK([XZ], [xz], [
              xz_version=`$ac_path_XZ --version 2>&1 | cut -d' ' -f4 | $SED -n 1p`
              AS_IF([test -n "$xz_version"], [
                  AX_COMPARE_VERSION([$xz_version], [ge], [4.999.0], [
                      ac_cv_path_XZ="$ac_path_XZ"
                  ])
              ])
            ])
            AS_IF([test -z "$ac_cv_path_XZ"], [sage_spkg_install_xz=yes])
        ])
    else
        sage_spkg_install_xz=yes
    fi
])
