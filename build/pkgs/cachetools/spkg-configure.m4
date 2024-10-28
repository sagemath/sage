SAGE_SPKG_CONFIGURE([cachetools], [
    sage_spkg_install_cachetools=yes
  ], [dnl REQUIRED-CHECK
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_TOX])
    dnl cachetools is only needed when we cannot use system tox.
    AS_VAR_SET([SPKG_REQUIRE], [$sage_spkg_install_tox])
  ])
