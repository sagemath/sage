.. _installation-guide:

Welcome to the SageMath Installation Guide
==========================================

You can install SageMath either from a package manager, a pre-built binary
tarball or from its sources.

Installing SageMath from your distribution package manager is the preferred and
fastest solution (dependencies will be automatically taken care of and SageMath
will be using your system Python). It is the case at least for the following
GNU/Linux distributions: Debian version >= 9, Ubuntu version >= 18.04,
Arch Linux, and NixOS. If you are in this situation, see
:ref:`sec-GNU-Linux`.

If your operating system does not provide SageMath, you can also use a
pre-built binary. See the section :ref:`sec-installation-from-binaries`.

Or you could install the ``sage`` package from the `conda-forge
<https://conda-forge.org/>`_ project. See the section
:ref:`sec-installation-conda`.

By compiling SageMath from its sources you might be able to run a slightly more
up-to-date version. You can also modify it and contribute back to the project.
Compiling SageMath might take up to 4 hours on a recent computer. To build
SageMath from source, go to the section :ref:`sec-installation-from-sources`.

Note that there are other alternatives to use SageMath that completely avoid
installing it:

- the `Sage Debian Live USB key <https://sagedebianlive.metelu.net/>`_: a full
  featured USB key that contains a whole Linux distribution including SageMath.
  This might be an option if you fail installing SageMath on your operating
  system.

- `CoCalc <https://cocalc.com/>`_: an online service that provides SageMath and
  many other tools.

- `Sage Cell Server <https://sagecell.sagemath.org/>`_: an online service for
  elementary SageMath computations.

- `Docker images <https://hub.docker.com/r/sagemath/sagemath/>`_: SageMath in a
  container for more experienced users.

The rest of this document describes how to install SageMath from pre-built
binaries and from sources.

.. toctree::
   :maxdepth: 2

   linux
   binary
   conda
   source
   launching
   troubles

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License`__.

__ http://creativecommons.org/licenses/by-sa/3.0/
