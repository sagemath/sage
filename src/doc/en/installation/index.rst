.. _installation-guide:

==================================
Welcome to Sage Installation Guide
==================================

Do you know that you can use Sage without installation? There are cloud
services that offer free access to Sage:

- `Sage Binder repo <https://github.com/sagemath/sage-binder-env>`_ provides a
  Binder badge to launch JupyterLab environment with Sage.

- `Sage Cell Server <https://sagecell.sagemath.org/>`_ is a free online service for
  quick computations with Sage.

- `CoCalc <https://cocalc.com/>`_ is an online commercial computing service. It
  offers resource-limited free access to Sage.

If none of these sounds interesting to you, you are welcome! Please read on.

This installation guide offers comprehensive information on installing Sage,
of the version specified on the left, below the logo, or of any recent
versions. However, you may need to check additional information and details
regarding installation on supported platforms available in the section
"Availability and installation help" of the `Release Tour
<https://github.com/sagemath/sage/releases>`_ for the Sage version you plan to
install.

There are alternative ways you can install Sage depending on how you plan to use it
and depending on your platform. So start in the relevant section below.

.. _installation-guide-windows:

Windows
-------

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

Start Ubuntu from the Start menu. You will see the command line prompt ``$``.

**I want to do Sage development:** Then follow the guide for Linux below.

**I don't want to do Sage development:** Then follow instructions in
:ref:`sec-installation-conda` and :ref:`sec-installation-conda-binary`. This is
the easiest way of using Sage in the terminal. Then see :ref:`sec-post-installation` for
recommended next steps, in particular for setting up the Jupyter notebook,
which is convenient if you want to use graphics.

.. _installation-guide-linux:

Linux
-----

**I want to do Sage development:**

- Obtain the Sage sources via ``git`` as described in `the Sage Developer Guide
  <https://doc.sagemath.org/html/en/developer/walkthrough.html#chapter-walkthrough>`_.

- Then build Sage from source, choosing ``develop`` git branch, as described in
  section :ref:`sec-installation-from-sources`.

- Alternatively, follow the instructions in section
  :ref:`sec-installation-conda-develop`;
  these describe a method that gets all required
  packages, including Python packages, from conda-forge.

**I don't want to do Sage development:**

- **I have root access (can use sudo command):** Then the easiest way to install Sage is
  through a Linux distribution that provides it as a package.  Some
  Linux distributions have up-to-date versions of Sage,
  see `repology.org: sagemath
  <https://repology.org/project/sagemath/versions>`_ for an
  overview. See :ref:`sec-installation-linux` for additional information.

  You may need to consider upgrading your Linux distribution because a recent
  version of Sage is only available on a newer version of the distribution.

  If there is no up-to-date version for your Linux distribution, you may
  consider `conda`, as explained below.

- **I have no root access or I am on an older distribution:** Install Sage from
  the `conda-forge <https://conda-forge.org/>`_ project, as described in section
  :ref:`sec-installation-conda`.

- Even if you don't want to do Sage development, you may attempt to build Sage
  from source as described in section :ref:`sec-installation-from-sources`
  choosing ``master`` git branch of the source code.

.. _installation-guide-macos:

macOS
-----

**I want to do Sage development:** Then follow the guide for Linux above.

**I don't want to do Sage development:**

- Install the binary build of Sage from the 3-manifolds project. See section
  :ref:`sec-installation-mac`.

- Alternatively, install Sage from
  the `conda-forge <https://conda-forge.org/>`_ project, as described in section
  :ref:`sec-installation-conda`.

- Even if you don't want to do Sage development, you may attempt to build Sage
  from source as described in section :ref:`sec-installation-from-sources`
  choosing ``master`` git branch of the source code and using Homebrew package
  manager.

Installation
------------

.. toctree::
   :maxdepth: 2

   source
   binary
   conda
   pypi
   launching
   troubles

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License`__.

__ http://creativecommons.org/licenses/by-sa/3.0/
