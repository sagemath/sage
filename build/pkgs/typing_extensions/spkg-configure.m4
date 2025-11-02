SAGE_SPKG_CONFIGURE([typing_extensions], [dnl
  dnl Although typing_extensions started as a backport package for Python <3.11,
  dnl it continues to provide additional typing features beyond what's in the stdlib
  dnl and is required as a runtime dependency by packages like beautifulsoup4.
  dnl Therefore, we always check for it regardless of Python version.
  SAGE_PYTHON_PACKAGE_CHECK([typing_extensions])
])
