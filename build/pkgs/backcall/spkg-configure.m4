SAGE_SPKG_CONFIGURE([backcall], [
  SAGE_PYTHON_PACKAGE_CHECK([backcall])
],[
  # required-check phase; skip this package if ipython
  # from the system is used (it's the only reverse dep).
  AC_REQUIRE([SAGE_SPKG_CONFIGURE_IPYTHON])
  AS_IF([test "${sage_spkg_install_ipython}" = "no"],[
    sage_require_backcall=no
  ])
])
