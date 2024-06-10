SAGE_SPKG_CONFIGURE([appdirs], [
    sage_spkg_install_appdirs=yes
  ], [dnl REQUIRED-CHECK
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_VIRTUALENV])
    dnl only needed as a dependency of virtualenv.
    AS_VAR_SET([SPKG_REQUIRE], [$sage_spkg_install_virtualenv])
  ])
