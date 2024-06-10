SAGE_SPKG_CONFIGURE([pdf2svg], [
   AC_PATH_PROG([PDF2SVG], [pdf2svg])
   AS_IF([test -z "$ac_cv_path_PDF2SVG"], [sage_spkg_install_pdf2svg=yes])
])

