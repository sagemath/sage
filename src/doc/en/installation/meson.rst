.. _build-source-meson:

================================
Building from source using Meson
================================

This is a short guide on how to build the Sage from source using Meson.

Walkthrough
===========

.. note::

    If you have previously build Sage in-place, you first have to delete the
    already compiled files, e.g. with ``shopt -s globstar`` followed by
    ``rm src/**/*.so`` or ``for f in src/**/*.so ; do mv "$f" "$f.old"; done``.
    Moreover, remove the old generated files with
    ``find src/sage/ext/interpreters -type f ! -name 'meson.build' -delete``.
    Also uninstall the 'old' sage packages with ``pip uninstall sage-conf sage-setup sagemath-standard``.

Sage relies on a number of external libraries, which have to be installed
before building. The easiest way to install them is to use Conda.
Alternatively, you can install them via your system package manager.
Both methods are described below.

Using Conda
~~~~~~~~~~~

- If you haven't already, download and install Miniforge for your platform
  following the `Miniforge instructions <https://github.com/conda-forge/miniforge?tab=readme-ov-file#install>`_.
  Other Conda distributions like Miniconda should work as well, but
  may require additional configuration (see :ref:`sec-installation-conda`).

- Create and activate a new conda environment with the dependencies of Sage
  and a few additional developer tools:

    .. tab:: Linux

        .. code-block:: shell

            $ mamba env create --file environment-3.12-linux.yml --name sage-dev
            $ mamba activate sage-dev

    .. tab:: macOS

        .. code-block:: shell

            $ mamba env create --file environment-3.12-macos.yml --name sage-dev
            $ mamba activate sage-dev

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

            $ mamba env create --file environment-3.12-win.yml --name sage-dev
            $ conda activate sage-dev
            $ set LIB=%CONDA_PREFIX%\Library\lib;%LIB%
            $ set INCLUDE=%CONDA_PREFIX%\Library\include;%INCLUDE%

        Windows support is experimental and not fully working yet.
        In fact, the Sage prompt is not working at all, but you can use the Python
        prompt to run certain commands. For example, the following should work
        after building Sage:

        .. code-block:: python

            >>> from sage.rings.integer import Integer
            >>> Integer(5)
            5
            >>> Integer(5) + 2.0
            7.0

    A different Python version can be selected by replacing ``3.12`` with the
    desired version.

- To compile and install Sage in editable install, just use:
  .. code-block:: shell-session

      $ pip install --no-build-isolation --editable .

  This will install Sage in the current Conda environment.
  The ``--no-build-isolation`` flag is necessary to allow the build system
  to reuse the already installed build dependencies.

- You can then start Sage from the command line with ``./sage``
  or run the tests with ``./sage -t``.

Using system package manager
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also install the dependencies via your system package manager.
Note, however, that not all dependencies may be available for your system,
and that the versions may be outdated.
In this case, Meson will try to build certain dependencies from source,
or it will fail with an error message.
In this case, you can either install the missing dependencies manually,
or use Conda to install them.

Depending on your distribution, install the following packages:

.. tab:: Debian/Ubuntu

    Not yet supported.

   .. .. literalinclude:: debian.txt

.. tab:: Fedora

    At least Fedora 41 is required.

   .. literalinclude:: fedora.txt


.. tab:: Arch Linux

   .. literalinclude:: arch.txt

.. tab:: Void Linux

   .. literalinclude:: void.txt


In the case that you want to install some dependencies manually, set the
correct environment variables to point to the installed libraries:

.. code-block:: shell-session

    $ export C_INCLUDE_PATH=$C_INCLUDE_PATH:/your/path/to/include
    $ export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/your/path/to/include
    $ export LIBRARY_PATH=$LIBRARY_PATH:/your/path/to/lib

We also recommend to install the Python package manager
`uv <https://docs.astral.sh/uv/getting-started/installation/>`_.

To compile and install Sage in editable install, then just use:

.. code-block:: shell-session

    $ uv venv
    $ uv pip install \
        meson-python \
        "cypari2 >=2.2.1" \
        "cython >=3.0, != 3.0.3, != 3.1.0" \
        "cython >=3.0, != 3.0.3" \
        "gmpy2 ~=2.1.b999" \
        memory_allocator \
        "numpy >=1.25" \
        jinja2 \
        setuptool
    $ uv sync --frozen --inexact --no-build-isolation

You can then start Sage from the command line with ``./sage``
or run the tests with ``./sage -t``.

Remarks
~~~~~~~

.. note::

    By using ``pip install --editable`` in the above steps, the Sage library
    is installed in editable mode. This means that when you only edit source
    files, there is no need to rebuild the library; it suffices to restart Sage.
    Note that this even works when you edit Cython files (they will be recompiled
    automatically), so you no longer need to manually compile after editing Cython
    files.

.. note::

    Note that ``make`` is not used at all, nor is ``configure``.
    This means that any Sage-the-distribution commands such as ``sage -i``
    will not work.

.. note::

    By default, Meson will automatically determine the number of jobs to
    run in parallel based on the number of CPU available. This can be adjusted
    by passing ``--config-settings=compile-args=-jN`` to ``pip install``.

    ``--verbose`` can be passed to ``pip install``, then the meson commands
    internally used by pip will be printed out.

.. note::

    To build the documentation, use::

        $ pip install --no-build-isolation -v -v --editable ./pkgs/sage-docbuild
        $ sage --docbuild all html

.. note::

  You can update the conda lock files by running ``tools/update-conda.py``.
  In order to update the conda environment afterwards use::

    $ mamba env update --file environment-3.12-linux.yml --name sage-dev


Background information
======================

Under the hood, pip invokes meson to configure and build the project.
We can also use meson directly as follows.

To configure the project, we need to run the following command:

.. code-block:: shell-session

    $ meson setup builddir

This will create a build directory ``builddir`` that will hold the build
artifacts.

If pip is used as above with ``--editable``, ``builddir`` is set to be
``build/cp[Python major version][Python minor version]``, such as
``build/cp311``.

To compile the project, run the following command:

.. code-block:: shell-session

    $ meson compile -C builddir

Installing is done with the following command:

.. code-block:: shell-session

    $ meson install -C builddir

This will install the project to currently active Python environment,
or to the system Python environment if no environment is active.
When editable install is used, it is not necessary to reinstall after each
compilation.

.. note::

    If you want to install the project to a different directory, you can specify
    the ``--prefix`` option when running the ``meson setup`` command:

    .. code-block:: shell-session

        $ meson setup builddir --prefix=/desired/install/path -Dpython.install_env=prefix

    This will then install in the directory specified by ``--prefix``,
    in particular the root folder will be be installed to
    ``/desired/install/path/lib/python3.12/site-packages/sage``.
    Usually, this directory is not on your Python path, so you have to use:

    .. code-block:: shell-session

        $ PYTHONPATH=/desired/install/path ./sage

Alternatively, we can still use pip to install:

.. code-block:: shell-session

    $ pip install --no-build-isolation --config-settings=builddir=builddir --editable .

.. tip::

    Package maintainers may want to specify further build options or need
    to install to a different directory than the install prefix.
    Both are supported naturally by Meson:

    .. code-block:: shell-session

        $ meson setup builddir --prefix=/usr --libdir=... -Dcpp_args=...
        $ meson compile -C builddir
        $ DESTDIR=/path/to/staging/root meson install -C builddir

    With the `default <https://mesonbuild.com/Running-Meson.html#installing>`_ prefix
    being ``/usr/local``, it may then install to
    ``$DESTDIR/usr/local/lib/python3.12/site-packages/sage``.

    See `Meson's quick guide <https://mesonbuild.com/Quick-guide.html#using-meson-as-a-distro-packager>`_
    and `Meson's install guide <https://mesonbuild.com/Installing.html#destdir-support>`_
    for more information.

Miscellaneous tips
==================

The environment variable ``MESONPY_EDITABLE_VERBOSE=1`` can be set while running ``./sage``,
so that when Cython files are recompiled a message is printed out.
See `<https://mesonbuild.com/meson-python/how-to-guides/editable-installs.html#verbose-mode>`_.

If a new ``.pyx`` file is added, it need to be added to ``meson.build`` file in
the containing directory.

Unlike the ``make``-based build system which relies on header comments
``# distutils: language = c++`` to determine whether C++ should be used,
Meson-based build system requires specifying
``override_options: ['cython_language=cpp']`` in the ``meson.build`` file.
Similarly, dependencies need to be specified by ``dependencies: [...]``.
