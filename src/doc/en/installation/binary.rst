.. _sec-installation-from-binaries:

Install from Prebuilt Binaries
==============================

.. _sec-installation-linux:

Linux
-----

Sage is available from various distributions and can be installed
by package managers.

As of Sage 10.2, we can recommend the following distributions, which
provide well-maintained and up-to-date Sage packages:
`Arch Linux <https://archlinux.org/>`_
and `Void Linux <https://voidlinux.org/>`_.

Gentoo users might want to give a try to
`sage-on-gentoo <https://github.com/cschwan/sage-on-gentoo>`_.

**Do not install a version of Sage older than 9.5.**
If you are on an older version of your distribution and a recent
version of Sage is only available on a newer version of the
distribution, consider upgrading your distribution.

See `the _sagemath dummy package <../reference/spkg/_sagemath.html>`_
for the names of packages that provide a standard installation of
Sage, including documentation and Jupyter.  See also `repology.org: sagemath
<https://repology.org/project/sagemath/versions>`_ for information
about versions of Sage packages in various distributions.

The `GitHub wiki page Distribution <https://github.com/sagemath/sage/wiki/Distribution>`_
collects information regarding packaging and distribution of Sage.

.. _sec-installation-mac:

macOS
-----

macOS binaries are available from `the 3-manifolds project <https://github.com/3-manifolds/Sage_macOS/releases/>`_.
It is a signed and notarized app, which works for macOS 10.12 and newer. It is
completely self-contained and provides the standard Sage distribution together
with many optional packages. Additional optional Python packages can be
installed with the ``%pip`` magic command and will go into your ``~/.sage``
directory.

Sage used to provide prebuilt binaries for macOS on its mirrors.
This has been discontinued, and the old binaries that are still available
there are no longer supported.

Windows
-------

Sage used to provide prebuilt binaries for Windows based on Cygwin.
This has been discontinued, and the old binaries that can be found
are no longer supported. Use Windows Subsystem for Linux instead.
