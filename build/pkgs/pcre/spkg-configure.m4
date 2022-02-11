SAGE_SPKG_CONFIGURE([pcre], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_BZIP2])
    AC_MSG_CHECKING([installing bzip2? ])
    if test x$sage_spkg_install_bzip2 = xyes; then
       AC_MSG_RESULT([yes; install pcre as well])
       sage_spkg_install_pcre=yes
    else
       AC_MSG_RESULT([no])
       dnl First try checking for pcre with pkg-config
       PKG_CHECK_MODULES([PCRE], [libpcre >= 8.39 libpcreposix libpcrecpp], [], [
           dnl Fallback to manually grubbing around for headers and libs
           AC_CHECK_HEADERS([pcre.h pcrecpp.h pcreposix.h], [
              AC_SEARCH_LIBS([regexec], [pcreposix], [], [sage_spkg_install_pcre=yes])],
           [sage_spkg_install_pcre=yes])
       ])
    fi
])

