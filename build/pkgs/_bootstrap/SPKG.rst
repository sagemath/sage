_bootstrap: Represents system packages required for running the top-level bootstrap script
==========================================================================================

Description
-----------

This optional script package represents the requirements (system packages)
that are needed in addition to those represented by the ``_prereq`` package
in order to run the top-level ``bootstrap`` script.

Namely, the following standard tools must be installed on your computer:

- **liblzma/xz**: XZ Utils, the free general-purpose data compression software
  with a high compression ratio.
- **python**: Python 3.12 or later from your system package manager or another
  external installation. Sage no longer builds its own Python interpreter; it
  uses this interpreter to create ``SAGE_VENV``. It needs to have the
  development headers and the following standard modules available:
  sqlite3, ctypes, math, hashlib, socket, ssl, ensurepip, zlib, setuptools

XZ Utils (liblzma) is available with most supported distributions (package
names such as ``xz-devel`` on Fedora, ``liblzma-dev`` on Debian/Ubuntu).
It can also be built from source from https://github.com/tukaani-project/xz.
After downloading and untarring the release archive, and changing to the
directory with the sources::

    $ ./configure --prefix=/usr/local && make && make install

will install XZ Utils in ``/usr/local``. Instead of ``/usr/local`` one may choose
another location, say ``/opt/foo``. A modern xz/liblzma installs a
``liblzma.pc`` file, and Sage's ``./configure`` checks for it using a
pkg-config macro. To detect xz/liblzma at a non-standard location, add it
to ``PKG_CONFIG_PATH``::

    $ export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/opt/foo/lib/pkgconfig

Python 3.12 (or later) with development headers is available with most
supported distributions (package names such as ``python3-devel`` on Fedora,
``python3-dev`` on Debian/Ubuntu).
It can also be built from source from https://www.python.org/downloads/.
After downloading and untarring the release archive, and changing to the
directory with the sources::

    $ ./configure --prefix=/usr/local --enable-optimizations --with-ensurepip=install && make && make install

will install Python in ``/usr/local``; this takes a few minutes on a
moderately fast machine. Instead of ``/usr/local`` one may choose another
location, say ``/opt/foo``, which then should be passed to Sage's
``./configure`` with ``--with-python=/opt/foo/bin/python3``.
The ``--enable-optimizations`` flag is optional but recommended for
performance. Make sure that the prerequisites for the required standard
modules (``sqlite3``, ``ctypes``, ``zlib``, ``ssl``, etc.) are installed
before building Python; otherwise these modules will be missing.

It is also possible to use a Python installed by `uv <https://docs.astral.sh/uv/>`_
and pass it to ``./configure`` with ``--with-python=...``.
