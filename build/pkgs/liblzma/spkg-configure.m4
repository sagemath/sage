SAGE_SPKG_CONFIGURE([liblzma], [
  SAGE_SPKG_DEPCHECK([xz], [
    dnl The library is actually installed by the xz spkg.
    AC_CHECK_LIB([lzma], [lzma_raw_decoder], [lzma_cv_liblzma=yes], [lzma_cv_liblzma=no])
    AC_CHECK_HEADER([lzma.h], [lzma_cv_lzma_h=yes], [lzma_cv_lzma_h=no])
    if test "$lzma_cv_liblzma" != "yes" -o "$lzma_cv_lzma_h" != "yes"; then
        sage_spkg_install_liblzma=yes
    fi
  ])
])
