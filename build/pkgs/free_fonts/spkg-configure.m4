SAGE_SPKG_CONFIGURE([free_fonts], [
    sage_spkg_install_free_fonts=yes
    m4_foreach([font],
               [FreeSerif.ttf,FreeSerif.otf],
        [
        AC_MSG_CHECKING([for ]font)
        AS_IF([kpsewhich ]font[ >& AS_MESSAGE_LOG_FD 2>&1], [
            AC_MSG_RESULT([yes])
            sage_spkg_install_free_fonts=no
            break
        ], [
            AC_MSG_RESULT([no])
        ])
    ])
])

