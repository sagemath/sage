.. _build-source-meson:

================================
Building from source using Meson
================================

This is a short guide on how to build the Sage from source using Meson.

Walkthrough
===========

Assume we're starting from a clean repo and a fully set up conda environment
(modify ``-linux`` according to your operating system):
        
.. tab:: Linux

    .. code-block:: shell

        $ mamba env create --file environment-3.11-linux.yml --name sage-dev
        $ conda activate sage-dev

.. tab:: macOS
    
    .. code-block:: shell

        $ mamba env create --file environment-3.11-macos.yml --name sage-dev
        $ conda activate sage-dev

.. tab:: Windows
    
    .. note::

        Windows support is very experimental and many features are not working
        yet.

    First you need to install the Microsoft Visual C++ compiler.
    You can download the 
    `Visual Studio Build Tools <https://aka.ms/vs/17/release/vs_BuildTools.exe>`_.
    Make sure to select "VC++ 2022 version xx.x build tools" and "Windows SDK".
    If you prefer, you can also run the following command to install the necessary
    components:
    
    .. code-block:: shell

        $ winget install Microsoft.VisualStudio.2022.BuildTools --force --override "--wait --passive --add Microsoft.VisualStudio.Component.VC.Tools.x86.x64 --add Microsoft.VisualStudio.Component.Windows11SDK.22621"
 
    Alternatively, you can use the compiler that comes bundled with Visual Studio.

    If you haven't already, install the latest version of Conda from
    `Miniforge <https://github.com/conda-forge/miniforge?tab=readme-ov-file#windows>`_.
    It is strongly recommended to choose the option to add Conda to the `PATH`
    during installation (because we will not use the Miniforge prompt).

    Open the "VS x64 Native Tools Command Prompt" (for 64bit) or 
    "Developer Command Prompt for VS2022 (or 2019)" (for 32bit).

    .. code-block:: shell

        $ mamba env create --file environment-3.11-win.yml --name sage-dev
        $ conda activate sage-dev
        $ set LIB=%CONDA_PREFIX%\Library\lib;%LIB%

    Windows support is experimental and not fully working yet.
    In fact, the Sage prompt is not working at all, but you can use the Python
    prompt to run certain commands. For example, the following should work:

    .. code-block:: python

        >>> from sage.rings.integer import Integer
        >>> Integer(5)
        5
        >>> Integer(5) + 2.0
        7.0


Alternatively, install all build requirements as described in section
:ref:`section-prereqs`. In the likely case that you have to install some
dependencies manually, set the correct environment variables to point
to the installed libraries:

.. CODE-BLOCK:: shell-session

    $ export C_INCLUDE_PATH=$C_INCLUDE_PATH:/your/path/to/include
    $ export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/your/path/to/include
    $ export LIBRARY_PATH=$LIBRARY_PATH:/your/path/to/lib

.. NOTE::

    If you have previously build Sage in-place, you first have to delete the
    already compiled files, e.g. with ``shopt -s globstar`` followed by 
    ``rm src/**/*.so`` or ``for f in src/**/*.so ; do mv "$f" "$f.old"; done``.
    Moreover, remove the old generated files with
    ``find src/sage/ext/interpreters -type f ! -name 'meson.build' -delete``. 
    Also uninstall the 'old' sage packages with ``pip uninstall sage-conf sage-setup sagemath-standard``.

To compile and install the project in editable install, just use:
    
.. CODE-BLOCK:: shell-session

    $ pip install --no-build-isolation --editable .

This will install Sage in the current Python environment. 
In a Conda environment, the ``--no-build-isolation`` flag is necessary to 
allow the build system to reuse the already installed build dependencies.
If you don't use Conda, you can omit this flag.

You can then start Sage from the command line with ``./sage`` 
or run the tests with ``./sage -t``.

.. NOTE::
    
    By using ``pip install --editable`` in the above steps, the Sage library 
    is installed in editable mode. This means that when you only edit source
    files, there is no need to rebuild the library; it suffices to restart Sage.
    Note that this even works when you edit Cython files (they will be recompiled
    automatically), so you no longer need to manually compile after editing Cython
    files.

.. NOTE::

    Note that ``make`` is not used at all, nor is ``configure``.
    This means that any Sage-the-distribution commands such as ``sage -i`` 
    will not work.

.. NOTE::

    By default, Meson will automatically determine the number of jobs to
    run in parallel based on the number of CPU available. This can be adjusted
    by passing ``--config-settings=compile-args=-jN`` to ``pip install``.

    ``--verbose`` can be passed to ``pip install``, then the meson commands
    internally used by pip will be printed out.

Background information
======================

Under the hood, pip invokes meson to configure and build the project.
We can also use meson directly as follows.

To configure the project, we need to run the following command:

.. CODE-BLOCK:: shell-session

    $ meson setup builddir --prefix=$PWD/build-install

This will create a build directory ``builddir`` that will hold the build artifacts.
The ``--prefix`` option specifies the directory where the Sage will be installed.

If pip is used as above, ``builddir`` is set to be
``build/cp[Python major version][Python minor version]``, such as ``build/cp311``.
``--prefix=`` can be left unspecified, when conda is used then meson will
install to the conda environment e.g. ``$HOME/miniforge3/envs/sage-dev/``.

To compile the project, run the following command:

.. CODE-BLOCK:: shell-session

    $ meson compile -C builddir

Installing is done with the following command:

.. CODE-BLOCK:: shell-session

    $ meson install -C builddir

This will then install in the directory specified by ``--prefix``, e.g.
``build-install/lib/python3.11/site-packages/sage``.
Usually, this directory is not on your Python path, so you have to use:

.. CODE-BLOCK:: shell-session

    $ PYTHONPATH=build-install/lib/python3.11/site-packages ./sage

When editable install is used, it is not necessary to reinstall after each compilation.

Alternatively, we can still use pip to install:

.. CODE-BLOCK:: shell-session

    $ pip install --no-build-isolation --config-settings=builddir=builddir --editable .

.. tip::

    Package maintainers may want to specify further build options or need
    to install to a different directory than the install prefix.
    Both are supported naturally by Meson:
    
    .. CODE-BLOCK:: shell-session

        $ meson setup builddir --prefix=/usr --libdir=... -Dcpp_args=...
        $ meson compile -C builddir
        $ DESTDIR=/path/to/staging/root meson install -C builddir
    
    See `Meson's quick guide <https://mesonbuild.com/Quick-guide.html#using-meson-as-a-distro-packager>`_
    and `Meson's install guide <https://mesonbuild.com/Installing.html#destdir-support>`_
    for more information.

Miscellaneous tips
==================

The environment variable ``MESONPY_EDITABLE_VERBOSE=1`` can be set while running ``./sage``,
so that when Cython files are recompiled a message is printed out.

If a new ``.pyx`` file is added, it need to be added to ``meson.build`` file in the
containing directory.

Unlike the ``make``-based build system which relies on header comments ``# distutils: language = c++``
to determine whether C++ should be used, Meson-based build system requires specifying
``override_options: ['cython_language=cpp']`` in the ``meson.build`` file.
Similarly, dependencies need to be specified by ``dependencies: [...]``.
