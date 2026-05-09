SAGE_SPKG_CONFIGURE([texlive_luatex], [
    sage_spkg_install_texlive_luatex=no
    AC_MSG_CHECKING([for luaotfload-main.lua])
    AS_IF([kpsewhich luaotfload-main.lua >& AS_MESSAGE_LOG_FD 2>&1], [
        AC_MSG_RESULT([yes])
    ], [
        AC_MSG_RESULT([no])
        sage_spkg_install_texlive_luatex=yes
    ])
])
