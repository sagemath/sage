SAGE_SPKG_CONFIGURE([pandoc], [
   AC_PATH_PROG([PANDOC], [pandoc])
   AS_IF([test -z "$ac_cv_path_PANDOC"], [sage_spkg_install_pandoc=yes])
])
