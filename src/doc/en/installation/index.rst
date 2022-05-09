.. _installation-guide:

=======================================
Welcome to the Sage Installation Guide!
=======================================

If you are reading this manual at https://doc.sagemath.org/, note that
it was built at the time the most recent stable release of SageMath
was made.

More up-to-date information and details regarding supported platforms
may have become available afterwards and can be found in the section
"Availability and installation help" of the
`release tour <https://wiki.sagemath.org/ReleaseTours>`_ for each
SageMath release.

**Where would you like to run SageMath?** Pick one of the following sections.

macOS
=====

- **Do you want to do SageMath development?**

  - **Yes, development:**

    Obtain the SageMath sources via ``git`` as described in `The Sage
    Developer's Guide
    <https://doc.sagemath.org/html/en/developer/walk_through.html#chapter-walkthrough>`_.

    - Then build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

    - Alternatively, follow the instructions in section
      :ref:`sec-installation-conda-develop`;
      these describe an experimental method that gets all required
      packages, including Python packages, from conda-forge.

  - **No development:**

    - Install the `binary build of SageMath <https://github.com/3-manifolds/Sage_macOS/releases>`_
      from the 3-manifolds project.  It is a signed and notarized app, which
      works for macOS 10.12 and newer. It is completely self-contained and
      provides the standard Sage distribution together with many optional
      packages. Additional optional Python packages can be installed with the
      ``%pip`` magic command and will go into your ``~/.sage`` directory.

    - Alternatively, install SageMath from the `conda-forge
      <https://conda-forge.org/>`_ project, as described in section
      :ref:`sec-installation-conda`.

    - Alternatively, build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

Windows
=======

- **Do you want to do SageMath development?**

  - **Yes, development:**

    Enable Windows Subsystem for Linux (WSL) by following the
    `official WSL setup guide
    <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_. Be
    sure to do the steps to install WSL2 and set it as default.
    Then go to the Microsoft Store and install Ubuntu (or another
    Linux distribution). Start Ubuntu from the start menu.

    Then follow the instructions for development on Linux below.

  - **No development:**

    - Enable Windows Subsystem for Linux (WSL) by following the
      `official WSL setup guide
      <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_. Be
      sure to do the steps to install WSL2 and set it as default.
      Then go to the Microsoft Store and install Ubuntu (or another
      Linux distribution). Start Ubuntu from the start menu.

      On the Linux running on WSL, you always have root access, so you
      can use any of the installation methods described below for
      Linux.

    - Alternatively, in particular if you cannot use WSL, install
      `Cygwin <https://cygwin.com/>`_ and then build SageMath from source
      as described in section :ref:`sec-installation-from-sources`.

Linux
=====

- **Do you want to do SageMath development?**

  - **Yes, development:**

    Obtain the SageMath sources via ``git`` as described in `The Sage
    Developer's Guide
    <https://doc.sagemath.org/html/en/developer/walk_through.html#chapter-walkthrough>`_.

    - Then build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

    - Alternatively, follow the instructions in section
      :ref:`sec-installation-conda-develop`;
      these describe an experimental method that gets all required
      packages, including Python packages, from conda-forge.

  - No development: **Do you have root access (sudo)?**

    - **Yes, root access:** Then the easiest way to install SageMath is
      through a Linux distribution that provides it as a package.  Most
      major Linux distributions have up-to-date versions of SageMath,
      see `repology.org: sagemath
      <https://repology.org/project/sagemath/versions>`_ for an
      overview.  See :ref:`sec-GNU-Linux` for additional information.

      If you are on an older version of your distribution and a recent
      version of SageMath is only available on a newer version of the
      distribution, consider upgrading your distribution.
      In particular, do not install a version of Sage older than 9.2.

    - **No root access, or on an older distribution** Install SageMath from
      the `conda-forge <https://conda-forge.org/>`_ project, as described in section
      :ref:`sec-installation-conda`.

    - Alternatively, build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

In the cloud
============

- `CoCalc <https://cocalc.com/>`_: an online service that provides SageMath and
  many other tools.

- On any system that allows you to bring your own Docker images to run in
  a container:  Use the `Docker image sagemath/sagemath <https://hub.docker.com/r/sagemath/sagemath/>`_.

- `Sage Cell Server <https://sagecell.sagemath.org/>`_: an online service for
  elementary SageMath computations.


More information:

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
