setuptools: Build system for Python packages
============================================

Description
-----------

setuptools is a collection of enhancements to the Python distutils (for
Python 2.6 and up) that allow you to more easily build and distribute
Python packages, especially ones that have dependencies on other
packages.

Website: http://pypi.python.org/pypi/setuptools/

License
-------

PSF or ZPL. i.e Python Software Foundation License or Zope Public
License


Upstream Contact
----------------

-  Phillip J. Eby (distutils-sig@python org)

Dependencies
------------

-  python


Build Instructions/Changes
--------------------------

The following patches are in the patches subdirectory. The patches are
applied during the build process.

-  pkg_resources.py.patch: silence warning about permissions.

-  easy_install_lock.patch: lock the easy_install.pth file to allow
   simultaneous installation
