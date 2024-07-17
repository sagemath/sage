SAGE_SPKG_CONFIGURE([colorama], [
    sage_spkg_install_colorama=yes
  ], [dnl REQUIRED-CHECK
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_TOX])
    dnl colorama is only needed when we cannot use system tox.
    AS_VAR_SET([SPKG_REQUIRE], [$sage_spkg_install_tox])
  ])
