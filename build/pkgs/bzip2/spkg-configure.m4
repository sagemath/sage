SAGE_SPKG_CONFIGURE([bzip2], [
    AC_CHECK_HEADER(bzlib.h, [
      AC_SEARCH_LIBS([BZ2_bzCompress], [bz2], [
        AC_PATH_PROG([bzip2_prog], [bzip2])
        AS_IF([test x$bzip2_prog = x], [sage_spkg_install_bzip2=yes])
      ], [sage_spkg_install_bzip2=yes])
    ], [sage_spkg_install_bzip2=yes])
])
