SAGE_SPKG_CONFIGURE([matplotlib], [
  SAGE_SPKG_DEPCHECK([freetype libpng qhull], [
    SAGE_PYTHON_PACKAGE_CHECK([matplotlib])
  ])
])
