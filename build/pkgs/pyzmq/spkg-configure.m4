SAGE_SPKG_CONFIGURE([pyzmq], [
  SAGE_SPKG_DEPCHECK([zeromq], [
    SAGE_PYTHON_PACKAGE_CHECK([pyzmq])
  ])
])
