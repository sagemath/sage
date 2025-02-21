SAGE_SPKG_CONFIGURE([debugpy], [
  SAGE_PYTHON_PACKAGE_CHECK([debugpy])
], [
  dnl REQUIRED-CHECK
  dnl
  dnl Skip debugpy if ipykernel from the system will be used.
  dnl This allows downstream packagers to treat debugpy (a
  dnl somewhat problematic package) as optional.
  AC_REQUIRE([SAGE_SPKG_CONFIGURE_IPYKERNEL])
  sage_require_debugpy=yes
  AS_VAR_IF([sage_spkg_install_ipykernel], [no], [
    sage_require_debugpy=no
  ])
])
