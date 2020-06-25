SAGE_SPKG_CONFIGURE([sympow], [
   AC_PATH_PROG([SYMPOW], [sympow])
   AS_IF([test -z "$ac_cv_path_SYMPOW"], [sage_spkg_install_sympow=yes])
])
