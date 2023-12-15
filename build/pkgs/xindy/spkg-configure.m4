SAGE_SPKG_CONFIGURE([xindy], [
    sage_spkg_install_xindy=no
    AC_PATH_PROG([XINDY], [xindy])
    AS_IF([test -z "$ac_cv_path_XINDY"], [sage_spkg_install_free_fonts=yes])
])
