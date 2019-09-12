SAGE_SPKG_CONFIGURE([bzip2], [
    AC_CHECK_HEADER(bzlib.h, [], [sage_spkg_install_bzip2=yes])
    AC_SEARCH_LIBS([BZ2_bzCompress], [bz2], [], [sage_spkg_install_bzip2=yes])
    AC_CHECK_PROG(bzip2, [break], [sage_spkg_install_bzip2=yes])
])
