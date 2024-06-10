SAGE_SPKG_CONFIGURE([filelock], [
    sage_spkg_install_filelock=yes
  ], [dnl REQUIRED-CHECK
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_VIRTUALENV])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_TOX])
    dnl only needed as a dependency of tox and virtualenv.
    AS_IF([test $sage_spkg_install_virtualenv = no -a $sage_spkg_install_tox = no],
          AS_VAR_SET([SPKG_REQUIRE], [no]))
  ])
