SAGE_SPKG_CONFIGURE([libpng], [
      dnl First try checking for libpng with pkg-config
      PKG_CHECK_MODULES([LIBPNG], [libpng >= 1.2], [], [
        sage_spkg_install_libpng=yes
        dnl Yes, we *could* fallback to manually grubbing around for headers and libs as follows:
        dnl   AC_CHECK_HEADERS([png.h], [break], [sage_spkg_install_libpng=yes])
        dnl   AC_SEARCH_LIBS([png_get_io_ptr], [png], [], [sage_spkg_install_libpng=yes])
        dnl But 'matplotlib' and 'sagelib' rely on pkg-config to locate libpng.
        dnl So we would have to tell them about the libpng that we found,
        dnl for example by creating a facade .pc file like we do for BLAS.
      ])
])
