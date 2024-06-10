SAGE_SPKG_CONFIGURE([imagemagick], [
   AC_PATH_PROG([CONVERT], [convert])
   AS_IF([test -z "$ac_cv_path_CONVERT"], [sage_spkg_install_imagemagick=yes])
])
