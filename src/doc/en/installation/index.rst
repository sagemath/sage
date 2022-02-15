.. _installation-guide:

Welcome to the SageMath Installation Guide
==========================================

If you are reading this manual at https://doc.sagemath.org/, note that
it was built at the time the most recent stable release of SageMath
was made.

More up-to-date information and details regarding supported platforms
may have become available afterwards and can be found in the section
"Availability and installation help" of the
`release tour <https://wiki.sagemath.org/ReleaseTours>`_ for each
SageMath release.

**Where would you like to run SageMath?**

**On Linux:**

- **Do you want to do SageMath development?**

  - **Yes, development:**

    Obtain the SageMath sources via ``git`` as described in `The Sage
    Developer's Guide
    <https://doc.sagemath.org/html/en/developer/walk_through.html#chapter-walkthrough>`_.

    - Then build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

    - Alternatively, follow the instructions in the SageMath wiki on
      `Conda for Sage Developers <https://wiki.sagemath.org/Conda>`_;
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

    - **No root access:** Install SageMath from the `conda-forge
      <https://conda-forge.org/>`_ project, as described in section
      :ref:`sec-installation-conda`.

    - Alternatively, build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

**On macOS:**

- **Do you want to do SageMath development?**

  - **Yes, development:**

    Obtain the SageMath sources via ``git`` as described in `The Sage
    Developer's Guide
    <https://doc.sagemath.org/html/en/developer/walk_through.html#chapter-walkthrough>`_.

    - Then build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

    - Alternatively, follow the instructions in the SageMath wiki on
      `Conda for Sage Developers <https://wiki.sagemath.org/Conda>`_;
      these describe an experimental method that gets all required
      packages, including Python packages, from conda-forge.

  - **No development:**

    - Install the `binary build of SageMath <https://github.com/3-manifolds/Sage_macOS/releases>`_
      from the 3-manifolds project.  It is a signed and notarized app.
      You can choose between two options: A 800 MB disk image with the standard Sage distribution,
      and a 1 GB disk image that adds many optional packages.

    - Alternatively, install SageMath from the `conda-forge
      <https://conda-forge.org/>`_ project, as described in section
      :ref:`sec-installation-conda`.

    - Alternatively, build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

**On Windows:**

- **Do you want to do SageMath development?**

  - **Yes, development:**

    Enable Windows Subsystem for Linux (WSL), selecting any Linux
    distribution, for example Ubuntu. Then follow the instructions for
    development on Linux.

  - **No development:**

    - Enable Windows Subsystem for Linux (WSL), selecting any Linux
      distribution, for example Ubuntu. On the Linux running on WSL,
      you always have root access, so you can use any of the
      installation methods described above for Linux.

    - Alternatively, in particular if you cannot use WSL, install
      `Cygwin <https://cygwin.com/>`_ and then build SageMath from source
      as described in section :ref:`sec-installation-from-sources`.

**In the cloud:**

- `CoCalc <https://cocalc.com/>`_: an online service that provides SageMath and
  many other tools.

- On any system that can allows you to bring your own Docker images to run in
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
