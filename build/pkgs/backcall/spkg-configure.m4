SAGE_SPKG_CONFIGURE([backcall], [
  SAGE_PYTHON_PACKAGE_CHECK([backcall])
],[
  # required-check phase; skip this package if ipython
  # from the system is used (it's the only reverse dep).
  SAGE_SPKG_DEPCHECK([ipython], [sage_require_backcall=no], [])
])
