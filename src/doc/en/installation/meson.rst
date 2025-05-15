.. _build-source-meson:

================================
Building from source using Meson
================================

This is a short guide on how to build the Sage from source using Meson.

Walkthrough
===========

Assume we're starting from a clean repo and a fully set up conda environment
(modify ``-linux`` according to your operating system):
        
.. CODE-BLOCK:: shell-session

    $ mamba env create --file environment-3.11-linux.yml --name sage-dev
    $ conda activate sage-dev

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
The ``--prefix`` option specifies the directory where the Sage will be installed,
and can be omitted when ``pip`` is used to install as explained below.

If pip is used as above with ``--editable``, ``builddir`` is set to be
``build/cp[Python major version][Python minor version]``, such as ``build/cp311``.

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

Alternatively, we can still use pip to install (which does not require specifying
``--prefix`` in advance and automatically works with conda environment):

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

If a new ``.pyx`` file is added, it need to be added to ``meson.build`` file in the
containing directory.

Unlike the ``make``-based build system which relies on header comments ``# distutils: language = c++``
to determine whether C++ should be used, Meson-based build system requires specifying
``override_options: ['cython_language=cpp']`` in the ``meson.build`` file.
Similarly, dependencies need to be specified by ``dependencies: [...]``.
