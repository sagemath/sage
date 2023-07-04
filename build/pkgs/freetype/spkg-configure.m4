SAGE_SPKG_CONFIGURE([freetype], [
    SAGE_SPKG_DEPCHECK([gcc libpng], [
      dnl freetype versions are libtool's ones, cf trac #30014
      PKG_CHECK_MODULES([FREETYPE], [freetype2 >= 20.0], [], [sage_spkg_install_freetype=yes])
    ])
], [], [], [
    if test x$sage_spkg_install_freetype = xyes; then
      AC_SUBST(SAGE_FREETYPE_PREFIX, ['$SAGE_LOCAL'])
    else
      AC_SUBST(SAGE_FREETYPE_PREFIX, [''])
    fi
])


