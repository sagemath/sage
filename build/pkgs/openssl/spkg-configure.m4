SAGE_SPKG_CONFIGURE([openssl], [
  AX_CHECK_OPENSSL([], [
      sage_spkg_install_openssl=yes
      AC_MSG_WARN([Because your system does not have a suitable OpenSSL library,
Sage will install a prerelease version of OpenSSL from the 3.0 alpha series.
The OpenSSL Project Team indicates that this prerelease version has been provided
for testing ONLY.  It should NOT be used for security critical purposes.
We strongly recommend to install OpenSSL using the system package manager and
to re-run configure.])
  ])
])
