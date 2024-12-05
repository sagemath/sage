.. _installation-guide:

==================================
Welcome to Sage Installation Guide
==================================

If you are reading this manual at https://doc.sagemath.org/, note that
it was built at the time the most recent stable release of SageMath
was made.

More up-to-date information and details regarding supported platforms
may have become available afterwards and can be found in the section
"Availability and installation help" of the
`release tour for each SageMath release <https://github.com/sagemath/sage/releases>`_.

**Where would you like to run SageMath?** Pick one of the following sections.

.. _installation-guide-macos:

macOS
=====

**Do you want to do SageMath development?**

-   **Yes, development:**

    Obtain the SageMath sources via ``git`` as described in `The Sage
    Developer's Guide
    <https://doc.sagemath.org/html/en/developer/walkthrough.html#chapter-walkthrough>`_.

    - Then build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

    - Alternatively, follow the instructions in section
      :ref:`sec-installation-conda-develop`;
      these describe an experimental method that gets all required
      packages, including Python packages, from conda-forge.

-   **No development:**

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

.. _installation-guide-windows:

Windows
=======

SageMath does not run natively on Windows. Using SageMath on this
platform requires a virtualization method. Hence, as a first step:

- Make sure that hardware-assisted virtualization is enabled in
  the EFI or BIOS of your system. If in doubt, refer to your
  system's documentation for instructions on how to do this.


**Are you familiar with or interested in learning about using the Linux command line?**

- **Yes, Linux command line:**

  Then we recommend the method described in section
  :ref:`installation-guide-windows-wsl` below.

  This method is also a good choice if you are not interested in using
  Visual Studio Code and if you want to have a persistent Sage
  installation for use in several projects.

- **No, prefer to use Visual Studio Code:**

  Then we recommend the method described in section
  :ref:`installation-guide-windows-devcontainers` below.

  `Visual Studio Code <https://code.visualstudio.com/>`_ is a popular
  and powerful code editor developed by Microsoft and made available
  free of charge. We recommmend it because of a key feature, the `Visual
  Studio Code Dev Containers extension
  <https://code.visualstudio.com/docs/devcontainers/containers>`_, which
  provides a very convenient way to install and use SageMath running in
  a Linux container.

  We recommend this installation method in particular to users who are not already
  familiar with Windows Subsystem for Linux. It is also suitable for machines that
  cannot run WSL 2 (because of hardware constraints or policy set by system
  administrators) but can run Docker.

.. _installation-guide-windows-wsl:

Installing SageMath in the Windows Subsystem for Linux
------------------------------------------------------

We recommend this installation method in particular to users who are already
familiar with Linux or Windows Subsystem for Linux. It is also a good choice if
you want to have a persistent Sage installation for use in several projects.

Enable `Windows Subsystem for Linux (WSL)
<https://learn.microsoft.com/en-us/windows/wsl/>`_ and install
Ubuntu as follows.

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

- If your computer has less than 10GB of RAM, `change the WSL settings
  <https://learn.microsoft.com/en-us/windows/wsl/wsl-config#main-wsl-settings>`_
  to make at least 5GB of RAM available to WSL.

Start Ubuntu from the Start menu.

**Do you want to do SageMath development?**

-   **Yes, development:**

    Then follow the :ref:`instructions for development on Linux <installation-guide-linux>` below.

-   **No development:**

    Type the following commands to install Sage from conda-forge.
    (The ``$`` represents the command
    line prompt, don't type it!) The second step will ask a few questions,
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

.. _installation-guide-windows-devcontainers:

Using Visual Studio Code with SageMath dev containers
-----------------------------------------------------

- `Download and install VS Code <https://code.visualstudio.com/>`_.

  .. note::

     Although the source code of VS Code is open source (MIT
     License), the Dev Containers extension is only known to work with the
     Visual Studio Code product, but not with fully open source builds such
     as `VSCodium <https://vscodium.com/>`_.

- In VS Code, if the `Visual Studio Code Dev Containers
  <https://code.visualstudio.com/docs/devcontainers/containers>`_ extension
  is not already installed, click the "Extension" icon on the left
  (or press :kbd:`Ctrl` + :kbd:`Shift` + :kbd:`X`) to open a list of
  extensions. Search for "Dev Containers" and install it.

  The extension will walk you through the installation of Docker
  Desktop and activating Windows Subsystem for Linux (WSL 2).

**Do you want to do SageMath development?**

-   **Yes, development:**

    - In VS Code, `clone
      <https://code.visualstudio.com/docs/editor/versioncontrol#_cloning-a-repository>`_
      the `SageMath git repository
      <https://github.com/sagemath/sage>`_. For example, use the VS Code
      Command Palette (:kbd:`Ctrl` + :kbd:`Shift` + :kbd:`P`), enter the
      command "Git: Clone", hit :kbd:`Enter`, paste the URL https://github.com/sagemath/sage
      and hit :kbd:`Enter` again.

    - Use the command palette (:kbd:`Ctrl` + :kbd:`Shift` + :kbd:`P`),
      enter the command "Dev Containers: Reopen Folder in Container" , and hit :kbd:`Enter`.

      By clicking on "show log", you can see what it does; it will take a while.

    - To use Sage in a terminal, `open a new terminal in VS Code
      <https://code.visualstudio.com/docs/terminal/basics>`_ and type ``./sage``.

-   **No development:**

    - Open the folder of a project where you would like to use Sage in VS Code, or
      just create a new folder for your work with Sage.

    - Create a subfolder named ``.devcontainer``. Download this file from the
      SageMath GitHub repository and place it in this subfolder:

      - https://github.com/sagemath/sage/blob/develop/.devcontainer/downstream-conda-forge-latest/devcontainer.json

    - Use the command palette (:kbd:`Ctrl` + :kbd:`Shift` + :kbd:`P`),
      enter the command "Dev Containers: Reopen Folder in Container" , and hit :kbd:`Enter`.

      By clicking on "show log", you can see what it does; it will take a while.

    - To use Sage in a terminal, `open a new terminal in VS Code
      <https://code.visualstudio.com/docs/terminal/basics>`_ and type ``sage``.

To use Sage in a Jupyter notebook, click the "Extension" icon on the left (or press
:kbd:`Ctrl` + :kbd:`Shift` + :kbd:`X`) to open a list of
extensions. Install the "Jupyter" extension. In the command
palette, enter "Create: New Jupyter Notebook", and hit
:kbd:`Enter`. Click "Select Kernel" on the right (or press :kbd:`Ctrl` +
:kbd:`Alt` + :kbd:`Enter`), select SageMath, and hit :kbd:`Enter`.


.. _installation-guide-linux:

Linux
=====

**Do you want to do SageMath development?**

-   **Yes, development:**

    Obtain the SageMath sources via ``git`` as described in `The Sage
    Developer's Guide
    <https://doc.sagemath.org/html/en/developer/walkthrough.html#chapter-walkthrough>`_.

    - Then build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

    - Alternatively, follow the instructions in section
      :ref:`sec-installation-conda-develop`;
      these describe an experimental method that gets all required
      packages, including Python packages, from conda-forge.

-   No development: **Do you have root access (sudo)?**

    - **Yes, root access:** Then the easiest way to install SageMath is
      through a Linux distribution that provides it as a package.  Some
      Linux distributions have up-to-date versions of SageMath,
      see `repology.org: sagemath
      <https://repology.org/project/sagemath/versions>`_ for an
      overview.  See :ref:`sec-GNU-Linux` for additional information.

      If you are on an older version of your distribution and a recent
      version of SageMath is only available on a newer version of the
      distribution, consider upgrading your distribution.
      In particular, do not install a version of Sage older than 9.5.

    - **No root access, or on an older distribution:** Install SageMath from
      the `conda-forge <https://conda-forge.org/>`_ project, as described in section
      :ref:`sec-installation-conda`.

    - Alternatively, build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

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


More information
================

.. toctree::
   :maxdepth: 2

   linux
   binary
   conda
   source
   meson
   launching
   troubles

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License`__.

__ http://creativecommons.org/licenses/by-sa/3.0/
