SAGE_SPKG_CONFIGURE([iconv], [
    AM_ICONV
    if test x"$am_cv_func_iconv" != xyes; then
       sage_spkg_install_iconv=yes
    fi
])

