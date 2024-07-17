SAGE_SPKG_CONFIGURE([xz], [
        AC_CACHE_CHECK([for xz >= 4.999.0], [ac_cv_path_XZ], [
            AC_PATH_PROGS_FEATURE_CHECK([XZ], [xz], [
              xz_version=`$ac_path_XZ --version 2>&1 | cut -d' ' -f4 | $SED -n 1p`
              AS_IF([test -n "$xz_version"], [
                  AX_COMPARE_VERSION([$xz_version], [ge], [4.999.0], [
                      ac_cv_path_XZ="$ac_path_XZ"
                      ac_path_XZ_found=:
                  ])
              ])
            ])
            AS_IF([test -z "$ac_cv_path_XZ"], [sage_spkg_install_xz=yes])
        ])
], [dnl REQUIRED-CHECK
        dnl Issue #30948: All dependencies on "xz" are merely build-time dependencies
        dnl on the xz binary (for unpacking the tarball in sage_bootstrap.uncompress.tar_file
        dnl - when sage-bootstrap-python is so old that it cannot do that by itself).
        dnl Packages that depend on actual xz or liblzma should depend on the liblzma spkg.
        AS_IF(["$SAGE_BOOTSTRAP_PYTHON" -c "import lzma" 2>& AS_MESSAGE_LOG_FD], [
             sage_require_xz=no
        ])
])
