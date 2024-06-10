SAGE_SPKG_CONFIGURE([curl], [
    # Only install curl (and libcurl) if we cannot find version 7.22 or
    # later, since that is what R needs.
    AC_CACHE_CHECK([for curl 7.22], [ac_cv_path_CURL], [
        AC_PATH_PROGS_FEATURE_CHECK([CURL], [curl], [
          ${ac_path_CURL}-config --checkfor 7.22 >/dev/null 2>/dev/null &&
                    ac_cv_path_CURL=${ac_path_CURL} &&
                    ac_path_CURL_found=:
        ])
    ])
    AS_IF([test -z "$ac_cv_path_CURL"], [sage_spkg_install_curl=yes])
    LIBCURL_CHECK_CONFIG(, 7.22, , [sage_spkg_install_curl=yes])
    # Anaconda on macOS provides a libcurl with @rpath. Check that
    # linking produces a binary than can be run. (#27941)
    AC_CACHE_CHECK([whether programs linking to libcurl can be executed],
      [sage_libcurl_cv_lib_curl_executable],
      [
           _libcurl_save_cppflags=$CPPFLAGS
           CPPFLAGS="$LIBCURL_CPPFLAGS $CPPFLAGS"
           _libcurl_save_libs=$LIBS
           LIBS="$LIBCURL $LIBS"

           AC_RUN_IFELSE([AC_LANG_PROGRAM([[#include <curl/curl.h>]],[[
             curl_easy_setopt(NULL,CURLOPT_URL,NULL);
           ]])], sage_libcurl_cv_lib_curl_executable=yes, sage_libcurl_cv_lib_curl_executable=no, [
              dnl cross compiling. link only
              AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <curl/curl.h>]],[[
               curl_easy_setopt(NULL,CURLOPT_URL,NULL);
             ]])], sage_libcurl_cv_lib_curl_executable=yes, sage_libcurl_cv_lib_curl_executable=no)]
           )
      ])
    AS_IF([test "$sage_libcurl_cv_lib_curl_executable" = "no"], [sage_spkg_install_curl=yes])
])
