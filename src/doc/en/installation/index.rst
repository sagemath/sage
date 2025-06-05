.. _installation-guide:

==================================
Welcome to Sage Installation Guide
==================================

This is the installation guide for SageMath, a free open-source mathematics software system.
It is designed to help you install SageMath on your computer.
For options to run SageMath in the cloud, see the section :ref:`sec-cloud` below.

If you are planning to do development on SageMath, please refer instead to the
`Sage Developer's Guide <../developer/walkthrough.html>`_ for instructions on
obtaining the source code and building SageMath.

.. tab:: Linux

  .. tab:: Conda

    Install SageMath from
    the `conda-forge <https://conda-forge.org/>`_ project, as described in section
    :ref:`sec-installation-conda`.

  .. tab:: Arch Linux/Manjaro

    Install SageMath by running the following command in a terminal::

      $ sudo pacman -S sagemath
  
  .. tab:: Void Linux
  
    Install SageMath by running the following command in a terminal::

      $ xbps-install -S sagemath

  .. tab:: Nix

    Install SageMath by running the following command in a terminal::

      $ nix-env -iA nixpkgs.sage

  .. tab:: Gentoo

    Install SageMath from the `sage-on-gentoo <https://github.com/cschwan/sage-on-gentoo>`_
    overlay by running the following command in a terminal::

      $ emerge --noreplace eselect-repository && eselect repository enable sage-on-gentoo && emerge --sync
      $ emerge -av sage

  .. tab:: Other Linux distributions

    Not all Linux distributions provide an up-to-date binary package of SageMath.
    **Do not install a version of Sage older than 9.5.**
    Instead we recommend to install SageMath via Conda, as described in the
    corresponding section.

    If you are on an older version of your distribution and a recent
    version of SageMath is only available on a newer version of the
    distribution, consider upgrading your distribution.

.. tab:: macOS 

  - Install the `binary build of SageMath <https://github.com/3-manifolds/Sage_macOS/releases>`_
    from the 3-manifolds project.  It is a signed and notarized app, which
    works for macOS 10.12 and newer. It is completely self-contained and
    provides the standard Sage distribution together with many optional
    packages. Additional optional Python packages can be installed with the
    ``%pip`` magic command and will go into your ``~/.sage`` directory.

  - Alternatively, install SageMath from the `conda-forge
    <https://conda-forge.org/>`_ project, as described in section
    :ref:`sec-installation-conda`.

.. tab:: Windows

    Enable `Windows Subsystem for Linux (WSL)
    <https://learn.microsoft.com/en-us/windows/wsl/>`_ and install
    Ubuntu as follows.

    - Make sure that hardware-assisted virtualization is enabled in
      the EFI or BIOS of your system. If in doubt, refer to your
      system's documentation for instructions on how to do this.

    - `Run the WSL install command as administrator.
      <https://learn.microsoft.com/en-us/windows/wsl/setup/environment#get-started>`_
      This will install Ubuntu Linux.

      Note that the basic instructions in the linked article apply to
      up-to-date installations of Windows 10 and 11, but there are
      also links to the procedures for older builds of Windows 10.

    - If you had installed WSL previously or installed it using
      different instructions, `verify that you are running WSL 2
      <https://learn.microsoft.com/en-us/windows/wsl/install#check-which-version-of-wsl-you-are-running>`_.

    - `Set up your Linux username and password.
      <https://learn.microsoft.com/en-us/windows/wsl/setup/environment#set-up-your-linux-username-and-password>`_
      Do not include any spaces in your username.

    - If your computer has less than 8GB of RAM, `change the WSL settings
      <https://learn.microsoft.com/en-us/windows/wsl/wsl-config#main-wsl-settings>`_
      to make at least 4GB of RAM available to WSL.

    Start Ubuntu from the Start menu, and type the following commands
    to install Sage from conda-forge. The second step will ask a few questions,
    and you may need to hit :kbd:`Enter` to confirm or type ``yes``
    and then hit :kbd:`Enter`.

    .. code-block:: shell

      $ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
      $ bash Miniforge3-$(uname)-$(uname -m).sh
      $ conda create -n sage sage python=3.11

    (If there are any installation failures, please report them to
    the conda-forge maintainers by opening a `GitHub Issue for
    conda-forge/sage-feedstock <https://github.com/conda-forge/sage-feedstock/issues>`_.)

    You can now start SageMath as follows:

    .. code-block:: shell

      $ conda activate sage
      $ sage

    This way of starting Sage gives you the most basic way of using
    Sage in the terminal. See :ref:`sec-launching` for recommended next steps,
    in particular for setting up the Jupyter notebook, which is required if
    you want to use graphics.

.. _sec-cloud:

In the cloud
============

- `Sage Binder repo <https://github.com/sagemath/sage-binder-env>`_ provides a
  Binder badge to launch JupyterLab environment with Sage.

- `Sage Cell Server <https://sagecell.sagemath.org/>`_ is a free online service for
  quick computations with Sage.

- `CoCalc <https://cocalc.com/>`_ is an online commercial service that provides Sage and
  many other tools.

- `Docker image sagemathinc/cocalc
  <https://hub.docker.com/r/sagemathinc/cocalc>`_ can be used on any system with Docker to run CoCalc locally.


More information:

.. toctree::
   :maxdepth: 2

   conda
   source
   meson
   launching
   troubles
