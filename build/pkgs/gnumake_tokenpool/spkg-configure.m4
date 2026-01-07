SAGE_SPKG_CONFIGURE([gnumake_tokenpool], [
  SAGE_PYTHON_PACKAGE_CHECK([gnumake_tokenpool])
], [dnl REQUIRED-CHECK
  dnl Do not install on macOS
  sage_require_gnumake_tokenpool=yes
  AS_CASE([$host], [*-*-darwin*], [sage_require_gnumake_tokenpool=no])
])
