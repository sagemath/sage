SAGE_SPKG_CONFIGURE([texlive], [
    sage_spkg_install_texlive=no
    AC_PATH_PROG([PDFLATEX], [pdflatex])
    AS_IF([test -z "$PDFLATEX"], [sage_spkg_install_texlive=yes])
    AC_PATH_PROG([LATEXMK], [latexmk])
    AS_IF([test -z "$LATEXMK"], [sage_spkg_install_texlive=yes])
    AC_PATH_PROG([DVIPNG], [dvipng])
    AS_IF([test -z "$DVIPNG"], [sage_spkg_install_texlive=yes])
    m4_foreach([latex_package],
               [fontspec,xunicode,xltxtra,amssymb,amsfonts,graphicx,mathrsfs,
                textcomp,tikz,tikz-qtree,iftex,tkz-berge,tkz-graph,xy,babel,
                subfigure,hyperref,hypcap,xr,tgtermes,fncychap],
        [
        AC_MSG_CHECKING([for latex package ]latex_package)
        AS_IF([kpsewhich ]latex_package[.sty >& AS_MESSAGE_LOG_FD 2>&1], [
            AC_MSG_RESULT([yes])
        ], [
            AC_MSG_RESULT([no])
            sage_spkg_install_texlive=yes
        ])
    ])
])
