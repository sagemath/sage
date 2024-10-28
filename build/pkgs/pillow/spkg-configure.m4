SAGE_SPKG_CONFIGURE([pillow], [
  SAGE_SPKG_DEPCHECK([bzip2 freetype libpng zlib], [
    SAGE_PYTHON_PACKAGE_CHECK([pillow])
  ])
])
