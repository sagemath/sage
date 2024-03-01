SAGE_SPKG_CONFIGURE([openssl], [
  AX_CHECK_OPENSSL([
    AC_MSG_CHECKING([whether OpenSSL >= 1.1.1, as required by PEP 644, and provides required APIs])
    AC_COMPILE_IFELSE(
        dnl Issue #32580: Need OpenSSL >= 1.1.1 for PEP 644
        dnl From https://www.openssl.org/docs/man3.0/man3/OPENSSL_VERSION_NUMBER.html:
        dnl If M  is the number from OPENSSL_VERSION_MAJOR
        dnl    NN is the number from OPENSSL_VERSION_MINOR
        dnl    PP is the number from OPENSSL_VERSION_PATCH
        dnl -> OPENSSL_VERSION_NUMBER is 0xMNN00PP0L
        dnl From https://www.openssl.org/docs/man1.1.1/man3/OPENSSL_VERSION_NUMBER.html
        dnl    FF is "fix"
        dnl    S  is "status" (f = release)
        dnl -> OPENSSL_VERSION_NUMBER is 0xMNNFFPPSL
        dnl
        dnl Issue #34273: Test program from ​https://github.com/python/cpython/blob/3.10/configure.ac#L5845
        [AC_LANG_PROGRAM([[
            #include <openssl/opensslv.h>
            #include <openssl/evp.h>
            #include <openssl/ssl.h>
            #if OPENSSL_VERSION_NUMBER < 0x10101000L
            #error OpenSSL >= 1.1.1 is required according to PEP 644
            #endif
            static void keylog_cb(const SSL *ssl, const char *line) {}
          ]], [[
            /* SSL APIs */
            SSL_CTX *ctx = SSL_CTX_new(TLS_client_method());
            SSL_CTX_set_keylog_callback(ctx, keylog_cb);
            SSL *ssl = SSL_new(ctx);
            X509_VERIFY_PARAM *param = SSL_get0_param(ssl);
            X509_VERIFY_PARAM_set1_host(param, "python.org", 0);
            SSL_free(ssl);
            SSL_CTX_free(ctx);
            /* hashlib APIs */
            OBJ_nid2sn(NID_md5);
            OBJ_nid2sn(NID_sha1);
            OBJ_nid2sn(NID_sha3_512);
            OBJ_nid2sn(NID_blake2b512);
            EVP_PBE_scrypt(NULL, 0, NULL, 0, 2, 8, 1, 0, NULL, 0);
          ]])
        ], [
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
  AC_REQUIRE([SAGE_SPKG_CONFIGURE_CURL])
  dnl openssl is a dependency only of python3 and curl; so if we use system python3 and curl,
  dnl we do not require it. (In particular, we do not need a specific version.)
  AS_IF([test x$sage_spkg_install_python3 = xno -a x$sage_spkg_install_curl = xno], [
    sage_require_openssl=no
  ])
])
