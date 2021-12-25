SAGE_SPKG_CONFIGURE([ffmpeg], [
   AC_PATH_PROG([FFMPEG], [ffmpeg])
   AS_IF([test -z "$ac_cv_path_FFMPEG"], [sage_spkg_install_ffmpeg=yes])
])

