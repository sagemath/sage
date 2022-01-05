SAGE_SPKG_CONFIGURE([virtualenv], [
    sage_spkg_install_virtualenv=yes
  ], [dnl REQUIRED-CHECK
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_TOX])
    dnl virtualenv is only needed when we cannot use system tox.
    AS_VAR_SET([SPKG_REQUIRE], [$sage_spkg_install_tox])
  ])
