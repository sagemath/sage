SAGE_SPKG_CONFIGURE([openssl], [
  AX_CHECK_OPENSSL([
    AC_MSG_CHECKING([whether OpenSSL >= 1.1.1, as required by PEP 644])
    AC_COMPILE_IFELSE(
        dnl Trac #32580: Need OpenSSL >= 1.1.1 for PEP 644
        dnl From https://www.openssl.org/docs/man3.0/man3/OPENSSL_VERSION_NUMBER.html:
        dnl If M  is the number from OPENSSL_VERSION_MAJOR
        dnl    NN is the number from OPENSSL_VERSION_MINOR
        dnl    PP is the number from OPENSSL_VERSION_PATCH
        dnl -> OPENSSL_VERSION_NUMBER is 0xMNN00PP0L
        dnl From https://www.openssl.org/docs/man1.1.1/man3/OPENSSL_VERSION_NUMBER.html
        dnl    FF is "fix"
        dnl    S  is "status" (f = release)
        dnl -> OPENSSL_VERSION_NUMBER is 0xMNNFFPPSL
        [AC_LANG_PROGRAM([[
            #include <openssl/ssl.h>
            #if OPENSSL_VERSION_NUMBER < 0x10101000L
            #  error OpenSSL >= 1.1.1 is required according to PEP 644
            #endif
        ]], [])], [
            AC_MSG_RESULT([yes])
            sage_spkg_install_openssl=no
        ], [
            AC_MSG_RESULT([no])
            sage_spkg_install_openssl=yes
        ])
  ], [dnl No openssl found
      sage_spkg_install_openssl=yes
  ])
], [dnl REQUIRED-CHECK
  AC_REQUIRE([SAGE_SPKG_CONFIGURE_PYTHON3])
  dnl openssl is a dependency only of python3; so if we use system python3,
  dnl we do not require it. (In particular, we do not need a specific version.)
  AS_IF([test x$sage_spkg_install_python3 = xno], [
    sage_require_openssl=no
  ])
])
