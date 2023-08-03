setuptools: Build system for Python packages
============================================

Description
-----------

setuptools is the classical build system for Python packages,
a collection of enhancements to the Python distutils.

This package represents version 63.x of ``setuptools``.
Sage installs this version to provide the build system
for non-PEP 517 packages. In particular, Sage uses it
for building ``numpy``, whose build system ``numpy.distutils`` 
is not compatible with newer versions of ``setuptools``,
see https://github.com/numpy/numpy/pull/22154

License
-------

MIT License

Upstream Contact
----------------

http://pypi.python.org/pypi/setuptools/

https://github.com/pypa/setuptools
