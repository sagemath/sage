SAGE_SPKG_CONFIGURE([ncurses], [
    dnl First try checking for ncurses with pkg-config
    PKG_CHECK_MODULES([NCURSES], [ncurses >= 6.0], [],
      [AC_CHECK_HEADERS([ncurses.h],
        [AC_SEARCH_LIBS([wresize], [ncurses tinfo], [],
                       [sage_spkg_install_ncurses=yes])],
        [sage_spkg_install_ncurses=yes])],
      [sage_spkg_install_ncurses=yes])
])
