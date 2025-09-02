SAGE_SPKG_CONFIGURE([pillow], [
  SAGE_SPKG_DEPCHECK([freetype libpng], [
    SAGE_PYTHON_PACKAGE_CHECK([pillow])
  ])
])
