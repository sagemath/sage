SAGE_SPKG_CONFIGURE([pillow], [
  SAGE_SPKG_DEPCHECK([freetype libpng zlib], [
    SAGE_PYTHON_PACKAGE_CHECK([pillow])
  ])
])
