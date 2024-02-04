SAGE_SPKG_CONFIGURE([appnope], [
  SAGE_PYTHON_PACKAGE_CHECK([appnope])
], [dnl REQUIRED-CHECK
  dnl Only required on macOS
  sage_require_appnope=no
  AS_CASE([$host], [*-*-darwin*], [sage_require_appnope=yes])
])
